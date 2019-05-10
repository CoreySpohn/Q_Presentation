#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu May  9 22:28:03 2019

@author: corey
"""
from big_ol_pile_of_manim_imports import *
PI = "3.14159265358979323846264338327950288419716939937510"
HAPPY = "Happy \pi Day!"

class RollAlongVector(Animation):
    CONFIG = {
        "rotation_vector" : OUT,
    }
    def __init__(self, mobject, vector, **kwargs):
        radius = mobject.get_width()/2
        radians = np.linalg.norm(vector)/radius
        last_alpha = 0
        digest_config(self, kwargs, locals())
        Animation.__init__(self, mobject, **kwargs)

    def update_mobject(self, alpha):
        d_alpha = alpha - self.last_alpha
        self.last_alpha = alpha
        self.mobject.rotate_in_place(
            -d_alpha*self.radians, 
            self.rotation_vector
        )
        self.mobject.shift(d_alpha*self.vector)


class PiDay(Scene):
    CYCLE = VGroup()
    def construct(self):
        self.divide_circle_into_parts()
        cyc = self.move_to_corner()
        self.rotate_wheels(cyc)

    def divide_circle_into_parts(self, parts=8):
        angle = 360/parts
        count = int(parts/2)
        circle = Circle(radius=2, color=WHITE)
        dot = Dot()
        start_line = Line(circle.get_left(),circle.get_right())
        self.play(ShowCreation(circle))
        self.play(ShowCreation(dot))
        self.play(Transform(dot, start_line))

        line_copy = start_line.copy()
        self.play(ShowCreation(start_line))

        for i in range(1, count):
            l_rotate = line_copy.rotate(angle=math.radians(angle), about_point=line_copy.get_center())
            to = Transform(start_line, l_rotate)
            self.play(to)
            self.add(l_rotate.copy())
#            self.dither()


        for i in self.get_mobjects():
            self.CYCLE.add(i)

    def move_to_corner(self):
        copy = self.CYCLE.copy().scale(0.5)
        copy.move_to(((LEFT*5+DOWN*2)))
        self.play(Transform(self.CYCLE,  copy))
        self.remove(self.CYCLE)
        return copy

    def rotate_wheels(self, group):
        cycle =  group.copy()

        bar = Line(cycle.get_center(), cycle.get_center()+UP*1.6)
        print(bar.get_end()+LEFT*0.5)
        print(bar.get_end()+RIGHT*0.5)
        handle = Mobject( Dot(cycle.get_center(),radius=0.15,color=RED),
                        Line(cycle.get_center(), cycle.get_center()+UP*1.6),
                        Line(bar.get_end()+LEFT*0.5,bar.get_end()+RIGHT*0.5))
        self.add(cycle)
        self.play(ShowCreation(handle))

        floor = Line(cycle.get_bottom()+LEFT*16, cycle.get_bottom()+RIGHT*16,color=BLACK)
        pi_floor = TexMobject(PI).scale(0.9).move_to((cycle.get_bottom()*1.063)+RIGHT*5.2)

        hpday = TexMobject(HAPPY, color=ORANGE).scale(0.9)
        rect = SurroundingRectangle(hpday)
        holder = Line(rect.get_bottom(), rect.get_bottom()*5)
        board = VGroup(hpday,rect, holder)

        pi = Randolph().scale(0.6)
        pi.move_to(handle.get_top()+UP*0.8)
        self.play(ShowCreation(pi))
        self.play(ShowCreation(board.move_to(pi.get_top()+RIGHT)))
    
        handle.add(pi)
        handle.add(board)
        self.add(floor)
        self.add(pi_floor)
        self.add(board)
        self.play(ShowCreation(pi_floor, run_time=4),RollAlongVector(cycle,floor.points[-1]-floor.points[0], rate_func=None, run_time=13),
            MaintainPositionRelativeTo(handle, cycle), rate_func=None)