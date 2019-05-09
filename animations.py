#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from big_ol_pile_of_manim_imports import *
class ControllabilityMatrix(Scene):
    def construct(self): 
        eq1=TextMobject(r"$Q_c = \left[\begin{matrix}0 & \frac{1}{l_{1}} & 0 & \frac{M g}{l_{1}^{2} m_{1}}\\\frac{1}{l_{1}} & 0 & \frac{M g}{l_{1}^{2} m_{1}} & 0\\0 & - \frac{1}{l_{2}} & 0 & - \frac{M g}{l_{2}^{2} m_{2}}\\- \frac{1}{l_{2}} & 0 & - \frac{M g}{l_{2}^{2} m_{2}} & 0\end{matrix}\right]$") 
        eq1.shift(1*UP)
        eq2=TextMobject(r"$rank(Q_c) = 4$")
        eq2.shift(1*DOWN)
        txt1=TextMobject("Therefore it is controllable")	
        txt1.shift(2*DOWN)
        self.play(Write(eq1))
        self.play(Write(eq2))
        self.play(FadeIn(txt1))
        self.wait(2)
        
class CartDemo(Scene):
    def construct(self):
        m1radius = 0.25
        m2radius = 0.25
        
        # Lengths, should also add angle
        l1 = 2
        l2 = 2
        
        # Angles (rad)
        t1 = 30 *np.pi/180
        t2 = 30 *np.pi/180
        
        # Unit vectors
        er1 = -np.sin(t1)*RIGHT+np.cos(t1)*UP
        et1 = -np.cos(t1)*RIGHT-np.sin(t1)*UP
        
        er2 = np.sin(t2)*RIGHT+np.cos(t2)*UP
        et2 = np.cos(t2)*RIGHT-np.sin(t2)*UP
        
        # Cart properties
        r_M_O = [0, -1] # Position of cart wrt the origin
        d_M = [1, 4] # Cart height, width
        
        # Wheel properties
        wheelradius = 0.4
        r_w1_O = [-2*wheelradius+r_M_O[0], -0.5*d_M[0]-wheelradius+r_M_O[1]]
        r_w2_O = [2*wheelradius+r_M_O[0], -0.5*d_M[0]-wheelradius+r_M_O[1]]
        
        # Ground properties
        r_G = [0,-0.5*d_M[0]-2*wheelradius-r_M_O[1]]
        
        # Position of hinges
        r_C1_O = [r_M_O[0]-0.25*d_M[1], r_M_O[1]+0.5*d_M[0], 0]
        r_C2_O = [r_M_O[0]+0.25*d_M[1], r_M_O[1]+0.5*d_M[0], 0]
        
        # Position of masses relative to their hinges
        r_m1_C1 = l1*er1
        r_m2_C2 = l2*er2
        
        # Position of masses relative to origin
        r_m1_O = [r_m1_C1[0]+r_C1_O[0], r_m1_C1[1]+r_C1_O[1], 0]
        r_m2_O = [r_m2_C2[0]+r_C2_O[0], r_m2_C2[1]+r_C2_O[1], 0]
        
        # Colors
        NEON_GREEN = '#39ff14'
        BLUE = '#4DA4B5'
        PURPLE = '#B696D2'
        MAROON = '#761648'
        GREY_1 = '#444444'
        
        # Create and place shapes
        M = Rectangle(height=d_M[0], width=d_M[1], fill_color=GREY_1, fill_opacity=1, stroke_color=NEON_GREEN)
        M.shift(r_M_O[0]*RIGHT+r_M_O[1]*UP)
        
        w1 = Circle(radius=wheelradius, fill_color=GREY_1, fill_opacity=1, stroke_color=NEON_GREEN)
        w1.shift(r_w1_O[0]*RIGHT+r_w1_O[1]*UP)
        
        w2 = Circle(radius=wheelradius, fill_color=GREY_1, fill_opacity=1, stroke_color=NEON_GREEN)
        w2.shift(r_w2_O[0]*RIGHT+r_w2_O[1]*UP)
        
        m1 = Circle(radius=m1radius, fill_color=GREY_1, fill_opacity=1, stroke_color=NEON_GREEN)
        m1.shift(r_m1_O[0]*RIGHT+r_m1_O[1]*UP)
        
        m2 = Circle(radius=m2radius, fill_color=GREY_1, fill_opacity=1, stroke_color=NEON_GREEN)
        m2.shift(r_m2_O[0]*RIGHT+r_m2_O[1]*UP)
          
        rod1 = Line(np.array(r_C1_O), np.array(r_m1_O), stroke_color=GREY_1)
        rod2 = Line(np.array(r_C2_O), np.array(r_m2_O), stroke_color=GREY_1)
        
        
        # Create groups to move together
        cart = VGroup(M, w1, w2)
        p1 = VGroup(m1, rod1)
        p2 = VGroup(m2, rod2)
        # Add elements to scene
        self.add(cart, p1, p2)
        self.bring_to_back(rod1, rod2)
        self.always_continually_update = True
        
        movement = 0.2*DOWN
        self.play(ApplyMethod(cart.move_to, self.relVec(cart, movement)), ApplyMethod(p1.move_to, self.relVec(p1, movement)), ApplyMethod(p2.move_to, self.relVec(p2, movement)))
#        self.play(FadeIn(m1), FadeIn(m2), FadeIn(rod1), FadeIn(rod2))
        self.bring_to_back(rod1, rod2)
        self.wait(2)
        
        
    def relVec(self, obj, pos):
        # generate the vector to move an object to a new position
        objpos = [obj.get_x(), obj.get_y(), obj.get_z()]
        return pos+objpos
            