import numpy as np
from pydub import AudioSegment
import shutil
import subprocess
import os
import _thread as thread
from time import sleep
import datetime

from manimlib.constants import FFMPEG_BIN
from manimlib.constants import STREAMING_IP
from manimlib.constants import STREAMING_PORT
from manimlib.constants import STREAMING_PROTOCOL
from manimlib.constants import VIDEO_DIR
from manimlib.utils.config_ops import digest_config
from manimlib.utils.file_ops import guarantee_existance
from manimlib.utils.file_ops import add_extension_if_not_present
from manimlib.utils.file_ops import get_sorted_integer_files
from manimlib.utils.sounds import get_full_sound_file_path


class SceneFileWriter(object):
    CONFIG = {
        "write_to_movie": False,
        # TODO, save_pngs is doing nothing
        "save_pngs": False,
        "png_mode": "RGBA",
        "save_last_frame": False,
        "movie_file_extension": ".mp4",
        "gif_file_extension": ".gif",
        "livestreaming": False,
        "to_twitch": False,
        "twitch_key": None,
        # Previous output_file_name
        # TODO, address this in extract_scene et. al.
        "file_name": None,
        "output_directory": None,
    }

    def __init__(self, scene, **kwargs):
        digest_config(self, kwargs)
        self.scene = scene
        self.stream_lock = False
        self.init_output_directories()
        self.init_audio()

    # Output directories and files
    def init_output_directories(self):
        output_directory = self.output_directory or self.get_default_output_directory()
        file_name = self.file_name or self.get_default_file_name()
        if self.save_last_frame:
            image_dir = guarantee_existance(os.path.join(
                VIDEO_DIR,
                output_directory,
                self.get_image_directory(),
            ))
            self.image_file_path = os.path.join(
                image_dir,
                add_extension_if_not_present(file_name, ".png")
            )
        if self.write_to_movie:
            movie_dir = guarantee_existance(os.path.join(
                VIDEO_DIR,
                output_directory,
                self.get_movie_directory(),
            ))
            self.movie_file_path = os.path.join(
                movie_dir,
                add_extension_if_not_present(
                    file_name, self.movie_file_extension
                )
            )
            self.gif_file_path = os.path.join(
                movie_dir,
                add_extension_if_not_present(
                    file_name, self.gif_file_extension
                )
            )
            self.partial_movie_directory = guarantee_existance(os.path.join(
                movie_dir,
                self.get_partial_movie_directory(),
                file_name,
            ))

    def get_default_output_directory(self):
        filename = os.path.basename(self.input_file_path)
        root, ext = os.path.splitext(filename)
        return root if root else ext[1:]

    def get_default_file_name(self):
        return self.scene.__class__.__name__

    def get_movie_directory(self):
        pixel_height = self.scene.camera.pixel_height
        frame_rate = self.scene.camera.frame_rate
        return "{}p{}".format(
            pixel_height, frame_rate
        )

    def get_image_directory(self):
        return "images"

    def get_partial_movie_directory(self):
        return "partial_movie_files"

    # Directory getters
    def get_image_file_path(self):
        return self.image_file_path

    def get_next_partial_movie_path(self):
        result = os.path.join(
            self.partial_movie_directory,
            "{:05}{}".format(
                self.scene.num_plays,
                self.movie_file_extension,
            )
        )
        return result

    def get_movie_file_path(self):
        return self.movie_file_path

    # Sound
    def init_audio(self):
        self.includes_sound = False

    def create_audio_segment(self):
        self.audio_segment = AudioSegment.silent()

    def add_audio_segment(self, new_segment,
                          time=None,
                          gain_to_background=None):
        if not self.includes_sound:
            self.includes_sound = True
            self.create_audio_segment()
        segment = self.audio_segment
        curr_end = segment.duration_seconds
        if time is None:
            time = curr_end
        if time < 0:
            raise Exception("Adding sound at timestamp < 0")

        new_end = time + new_segment.duration_seconds
        diff = new_end - curr_end
        if diff > 0:
            segment = segment.append(
                AudioSegment.silent(int(np.ceil(diff * 1000))),
                crossfade=0,
            )
        self.audio_segment = segment.overlay(
            new_segment,
            position=int(1000 * time),
            gain_during_overlay=gain_to_background,
        )

    def add_sound(self, sound_file, time=None, gain=None, **kwargs):
        file_path = get_full_sound_file_path(sound_file)
        new_segment = AudioSegment.from_file(file_path)
        if gain:
            new_segment = new_segment.apply_gain(gain)
        self.add_audio_segment(new_segment, time, **kwargs)

    # Writers
    def begin_animation(self, allow_write=False):
        if self.write_to_movie and allow_write:
            self.open_movie_pipe()
        if self.livestreaming:
            self.stream_lock = False

    def end_animation(self, allow_write=False):
        if self.write_to_movie and allow_write:
            self.close_movie_pipe()
        if self.livestreaming:
            self.stream_lock = True
            thread.start_new_thread(self.idle_stream, ())

    def write_frame(self, frame):
        if self.write_to_movie:
            self.writing_process.stdin.write(frame.tostring())

    def save_final_image(self, image):
        file_path = self.get_image_file_path()
        image.save(file_path)
        self.print_file_ready_message(file_path)

    def idle_stream(self):
        while self.stream_lock:
            a = datetime.datetime.now()
            self.update_frame()
            n_frames = 1
            frame = self.get_frame()
            self.add_frames(*[frame] * n_frames)
            b = datetime.datetime.now()
            time_diff = (b - a).total_seconds()
            frame_duration = 1 / self.scene.camera.frame_rate
            if time_diff < frame_duration:
                sleep(frame_duration - time_diff)

    def finish(self):
        if self.write_to_movie:
            if hasattr(self, "writing_process"):
                self.writing_process.terminate()
            self.combine_movie_files()
        if self.save_last_frame:
            self.scene.update_frame(ignore_skipping=True)
            self.save_final_image(self.scene.get_image())

    def open_movie_pipe(self):
        file_path = self.get_next_partial_movie_path()
        temp_file_path = os.path.splitext(file_path)[0] + '_temp' + self.movie_file_extension

        self.partial_movie_file_path = file_path
        self.temp_partial_movie_file_path = temp_file_path

        fps = self.scene.camera.frame_rate
        height = self.scene.camera.get_pixel_height()
        width = self.scene.camera.get_pixel_width()

        command = [
            FFMPEG_BIN,
            '-y',  # overwrite output file if it exists
            '-f', 'rawvideo',
            '-s', '%dx%d' % (width, height),  # size of one frame
            '-pix_fmt', 'rgba',
            '-r', str(fps),  # frames per second
            '-i', '-',  # The imput comes from a pipe
            '-c:v', 'h264_nvenc',
            '-an',  # Tells FFMPEG not to expect any audio
            '-loglevel', 'error',
        ]
        if self.movie_file_extension == ".mov":
            # This is if the background of the exported video
            # should be transparent.
            command += [
                '-vcodec', 'qtrle',
                # '-vcodec', 'png',
            ]
        else:
            command += [
                '-vcodec', 'libx264',
                '-pix_fmt', 'yuv420p',
            ]
        if self.livestreaming:
            if self.to_twitch:
                command += ['-f', 'flv']
                command += ['rtmp://live.twitch.tv/app/' + self.twitch_key]
            else:
                command += ['-f', 'mpegts']
                command += [STREAMING_PROTOCOL + '://' + STREAMING_IP + ':' + STREAMING_PORT]
        else:
            command += [temp_file_path]
        self.writing_process = subprocess.Popen(command, stdin=subprocess.PIPE)

    def close_movie_pipe(self):
        self.writing_process.stdin.close()
        self.writing_process.wait()
        if self.livestreaming:
            return True
        shutil.move(
            self.temp_partial_movie_file_path,
            self.partial_movie_file_path,
        )

    def combine_movie_files(self):
        # Manim renders the scene as many smaller movie files
        # which are then concatenated to a larger one.  The reason
        # for this is that sometimes video-editing is made easier when
        # one works with the broken up scene, which effectively has
        # cuts at all the places you might want.  But for viewing
        # the scene as a whole, one of course wants to see it as a
        # single piece.
        kwargs = {
            "remove_non_integer_files": True,
            "extension": self.movie_file_extension,
        }
        if self.scene.start_at_animation_number is not None:
            kwargs["min_index"] = self.scene.start_at_animation_number
        if self.scene.end_at_animation_number is not None:
            kwargs["max_index"] = self.scene.end_at_animation_number
        else:
            kwargs["remove_indices_greater_than"] = self.scene.num_plays - 1
        partial_movie_files = get_sorted_integer_files(
            self.partial_movie_directory,
            **kwargs
        )
        if len(partial_movie_files) == 0:
            print("No animations in this scene")
            return

        # Write a file partial_file_list.txt containing all
        # partial movie files
        file_list = os.path.join(
            self.partial_movie_directory,
            "partial_movie_file_list.txt"
        )
        with open(file_list, 'w') as fp:
            for pf_path in partial_movie_files:
                if os.name == 'nt':
                    pf_path = pf_path.replace('\\', '/')
                fp.write("file \'{}\'\n".format(pf_path))

        movie_file_path = self.get_movie_file_path()
        commands = [
            FFMPEG_BIN,
            '-y',  # overwrite output file if it exists
            '-f', 'concat',
            '-safe', '0',
            '-i', file_list,
            '-loglevel', 'error',
            
        ]
        if not self.save_as_gif:
            commands +=[
                '-c', 'copy',
                movie_file_path
            ]
        if self.save_as_gif:
            movie_file_path=self.gif_file_path
            commands +=[
                movie_file_path,
            ]
        if not self.includes_sound:
            commands.insert(-1, '-an')

        combine_process = subprocess.Popen(commands)
        combine_process.wait()

        if self.includes_sound:
            sound_file_path = movie_file_path.replace(
                self.movie_file_extension, ".wav"
            )
            # Makes sure sound file length will match video file
            self.add_audio_segment(AudioSegment.silent(0))
            self.audio_segment.export(
                sound_file_path,
                bitrate='312k',
            )
            temp_file_path = movie_file_path.replace(".", "_temp.")
            commands = [
                "ffmpeg",
                "-i", movie_file_path,
                "-i", sound_file_path,
                '-y',  # overwrite output file if it exists
                "-c:v", "copy",
                "-c:a", "aac",
                "-b:a", "320k",
                # select video stream from first file
                "-map", "0:v:0",
                # select audio stream from second file
                "-map", "1:a:0",
                '-loglevel', 'error',
                # "-shortest",
                temp_file_path,
            ]
            subprocess.call(commands)
            shutil.move(temp_file_path, movie_file_path)
            subprocess.call(["rm", sound_file_path])

        self.print_file_ready_message(movie_file_path)

    def print_file_ready_message(self, file_path):
        print("\nFile ready at {}\n".format(file_path))
