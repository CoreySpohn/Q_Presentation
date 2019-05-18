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

class Dynamics(Scene):
    def construct(self):
#        self.energy2L()
        self.l2EOM()
        
    def energy2L(self):
        # Kinetic Energy
        Teq = TexMobject(r"KE =")
        Teq1 = TexMobject(r"0.5 M \dot{x}^{2} + 0.5 m_{1} \left(l_{1}^{2} \dot{\theta}_{1}^{2} - 2 l_{1} \operatorname{cos}\left(\theta_{1}\right) \dot{\theta}_{1} \dot{x} + \dot{x}^{2}\right)")
        Teq2 = TexMobject(r"+ 0.5 m_{2} \left(l_{2}^{2} \dot{\theta}_{2}^{2} + 2 l_{2} \operatorname{cos}\left(\theta_{2}\right) \dot{\theta}_{2} \dot{x} + \dot{x}^{2}\right)")
        Teq.shift(2*UP+6*LEFT)
        Teq1.shift(2*UP)
        Teq2.shift(2*UP)
        Teq1.next_to(Teq, RIGHT)
        Teq2.next_to(Teq1, DOWN)
        Teq2.align_to(Teq1, LEFT)
        KEgroup = VGroup(Teq, Teq1, Teq2)
        
        # Potential Energy
        Ueq = TexMobject(r"PE = g \left(l_{1} m_{1} \operatorname{cos}\left(\theta_{1}\right) + l_{2} m_{2} \operatorname{cos}\left(\theta_{2}\right)\right)")
        Ueq.next_to(Teq2, DOWN)
        Ueq.align_to(Teq, LEFT)

        # Lagrangian symbolic
        Leq1_1 = TexMobject(r"\mathcal{L} =")
        Leq1_2 = TexMobject(r"KE - PE")
        Leq1_1.next_to(Ueq, DOWN)
        Leq1_1.align_to(Teq, RIGHT)
        Leq1_2.next_to(Leq1_1, RIGHT)
        Leq1group = VGroup(Leq1_1, Leq1_2)
        
        # Lagrangian full
        Leq2 = TexMobject(r"\mathcal{L} =")
        Leq2_1 = TexMobject(r"0.5 M \dot{x}^{2} + 0.5 m_{1} \left(l_{1}^{2} \dot{\theta}_{1}^{2} - 2 l_{1} \operatorname{cos}\left(\theta_{1}\right) \dot{\theta}_{1} \dot{x} + \dot{x}^{2}\right)")
        Leq2_2 = TexMobject(r"+ 0.5 m_{2} \left(l_{2}^{2} \dot{\theta}_{2}^{2} + 2 l_{2} \operatorname{cos}\left(\theta_{2}\right) \dot{\theta}_{2} \dot{x} + \dot{x}^{2}\right)")
        Leq2_3 = TexMobject(r"-g \left(l_{1} m_{1} \operatorname{cos}\left(\theta_{1}\right) + l_{2} m_{2} \operatorname{cos}\left(\theta_{2}\right)\right)")
        
        Leq2.move_to(Leq1group)
        Leq2.align_to(Leq1_1, LEFT)
        Leq2_1.next_to(Leq2, RIGHT)
        Leq2_2.next_to(Leq2_1, DOWN)
        Leq2_3.next_to(Leq2_2, DOWN)
        Leq2_2.align_to(Leq2_1, LEFT)
        Leq2_3.align_to(Leq2_1, LEFT)
        
        Lgroup = VGroup(Leq2, Leq2_1, Leq2_2, Leq2_3)
        
        # Kinetic and Potential Animations
        self.play(Write(KEgroup), Write(Ueq))
        self.wait(0.5)
        
        # Lagrangian Animations
        self.play(ApplyMethod(KEgroup.shift, 1*UP), ApplyMethod(Ueq.shift, 1*UP), FadeInFromDown(Leq1group))
        self.wait(0.5)
        self.play(ReplacementTransform(Leq1group, Lgroup))
        self.wait(0.5)
        
        # Focus on Lagrangian
        self.play(FadeOut(KEgroup), FadeOut(Ueq), ApplyMethod(Lgroup.center))
        self.remove(Lgroup)
        
    def l2EOM(self):
        # Lagrange equation 
        Leq2 = TexMobject(r"\mathcal{L} =")
        Leq2_1 = TexMobject(r"0.5 M \dot{x}^{2} + 0.5 m_{1} \left(l_{1}^{2} \dot{\theta}_{1}^{2} - 2 l_{1} \operatorname{cos}\left(\theta_{1}\right) \dot{\theta}_{1} \dot{x} + \dot{x}^{2}\right)")
        Leq2_2 = TexMobject(r"+ 0.5 m_{2} \left(l_{2}^{2} \dot{\theta}_{2}^{2} + 2 l_{2} \operatorname{cos}\left(\theta_{2}\right) \dot{\theta}_{2} \dot{x} + \dot{x}^{2}\right)")
        Leq2_3 = TexMobject(r"-g \left(l_{1} m_{1} \operatorname{cos}\left(\theta_{1}\right) + l_{2} m_{2} \operatorname{cos}\left(\theta_{2}\right)\right)")
        
        Leq2_1.next_to(Leq2, RIGHT)
        Leq2_2.next_to(Leq2_1, DOWN)
        Leq2_3.next_to(Leq2_2, DOWN)
        Leq2_2.align_to(Leq2_1, LEFT)
        Leq2_3.align_to(Leq2_1, LEFT)
        
        Lgroup = VGroup(Leq2, Leq2_1, Leq2_2, Leq2_3)
        Lgroup.center()
        self.play(Lgroup.to_edge, UP)
        
        # The Lagrange equations of motion
        Leq3 = TexMobject(r" \frac{d}{dt}\left(\frac{\partial \mathcal{L}}{\partial \dot{q}}\right) -\frac{\partial \mathcal{L}}{\partial q} = \sum \frac{\partial \vec{r}}{\partial q} \cdot \vec{F}")
        Leq3.next_to(Lgroup, DOWN)
        self.play(Write(Leq3))
        # Scale and move
        self.play(ApplyMethod(Lgroup.scale,0.75), ApplyMethod(Leq3.scale,0.75))
        self.play(ApplyMethod(Lgroup.to_corner,UP+LEFT))
        self.play(Leq3.next_to,Lgroup,RIGHT)
        
        # Generalized coordinates
        qall = TexMobject(r"q=(x, \theta_1, \theta_2)")
        qx = TexMobject(r"q=x :")
        qt1 = TexMobject(r"q=\theta_1 :")
        qt2 = TexMobject(r"q=\theta_2 :")
        
        qall.next_to(Lgroup,DOWN)
        qall.set_x(0)
        
        qx.to_edge(LEFT)
        qt1.to_edge(LEFT)
        qt2.to_edge(LEFT)
        
        qt1.shift(1.5*DOWN)
        qt2.shift(3*DOWN)
        
        # Equations for each q val
        xeq = TexMobject(r"M \ddot{x} + m_{1} \left(l_{1} \operatorname{sin}\left(\theta_{1}\right) \dot{\theta}_{1}^{2} - l_{1} \operatorname{cos}\left(\theta_{1}\right) \ddot{\theta}_{1} + \ddot{x}\right) + m_{2} \left(- l_{2} \operatorname{sin}\left(\theta_{2}\right) \dot{\theta}_{2}^{2} + l_{2} \operatorname{cos}\left(\theta_{2}\right) \ddot{\theta}_{2} + \ddot{x}\right) = u")
        xeq.scale(0.60)
        xeq.next_to(qx, RIGHT)

        t1eq = TexMobject(r"l_{1} m_{1} \left(- g \operatorname{sin}\left(\theta_{1}\right) + l_{1} \ddot{\theta}_{1} - \operatorname{cos}\left(\theta_{1}\right) \ddot{x}\right) = 0")
        t1eq.scale(0.60)
        t1eq.next_to(qt1, RIGHT)
        
        t2eq = TexMobject(r"l_{2} m_{2} \left(- g \operatorname{sin}\left(\theta_{2}\right) + l_{2} \ddot{\theta}_{2} + \operatorname{cos}\left(\theta_{2}\right) \ddot{x}\right) = 0")
        t2eq.scale(0.60)
        t2eq.next_to(qt2, RIGHT)
        
        self.play(Write(qall))
        self.play(TransformFromCopy(qall, qx), TransformFromCopy(qall, qt1), TransformFromCopy(qall, qt2))
        
        self.play(Write(xeq), Write(t1eq), Write(t2eq))
        self.wait(0.5)
        
        
        
        # Equations of motion
        xddLHS = TexMobject(r"\ddot{x} = ")
        xddRHS = TexMobject(r"\frac{0.5 g m_{1} \operatorname{sin}\left(2\theta_{1}\right) - 0.5 g m_{2} \operatorname{sin}\left(2\theta_{2}\right) - l_{1} m_{1} \operatorname{sin}\left(\theta_{1}\right) \dot{\theta}_{1}^{2} + l_{2} m_{2} \operatorname{sin}\left(\theta_{2}\right) \dot{\theta}_{2}^{2} + u}{M + m_{1} \operatorname{sin}^{2}\left(\theta_{1}\right) + m_{2} \operatorname{sin}^{2}\left(\theta_{2}\right)}")
        xddRHS.scale(0.45)
        xddLHS.next_to(qx, RIGHT)
        xddRHS.next_to(xddLHS, RIGHT)
        xdd = VGroup(xddLHS, xddRHS)
        xdd.set_x(0)
        
        t1ddLHS = TexMobject(r"\ddot{\theta_1} =")
        t1ddRHS = TexMobject(r"\frac{- 0.25 g m_{2} \left(- \operatorname{sin}\left(\theta_{1} - 2 \theta_{2}\right) + \operatorname{sin}\left(\theta_{1} + 2 \theta_{2}\right)\right) + g \left(M + m_{1} + m_{2} \operatorname{sin}^{2}\left(\theta_{2}\right)\right) \operatorname{sin}\left(\theta_{1}\right) + \left(- l_{1} m_{1} \operatorname{sin}\left(\theta_{1}\right) \dot{\theta}_{1}^{2} + l_{2} m_{2} \operatorname{sin}\left(\theta_{2}\right) \dot{\theta}_{2}^{2} + u\right) \operatorname{cos}\left(\theta_{1}\right)}{l_{1} \left(M + m_{1} \operatorname{sin}^{2}\left(\theta_{1}\right) + m_{2} \operatorname{sin}^{2}\left(\theta_{2}\right)\right)}")
        t1ddRHS.scale(0.45)
        t1ddLHS.next_to(qt1, RIGHT)
        t1ddRHS.next_to(t1ddLHS,RIGHT)
        t1dd = VGroup(t1ddLHS, t1ddRHS)
        t1dd.next_to(xdd, DOWN)
        
        t2ddLHS = TexMobject(r"\ddot{\theta_2} =")
        t2ddRHS = TexMobject(r"\frac{- 0.25 g m_{1} \left(\operatorname{sin}\left(2 \theta_{1} - \theta_{2}\right) + \operatorname{sin}\left(2 \theta_{1} + \theta_{2}\right)\right) + g \left(M + m_{1} \operatorname{sin}^{2}\left(\theta_{1}\right) + m_{2}\right) \operatorname{sin}\left(\theta_{2}\right) - \left(- l_{1} m_{1} \operatorname{sin}\left(\theta_{1}\right) \dot{\theta}_{1}^{2} + l_{2} m_{2} \operatorname{sin}\left(\theta_{2}\right) \dot{\theta}_{2}^{2} + u\right) \operatorname{cos}\left(\theta_{2}\right)}{l_{2} \left(M + m_{1} \operatorname{sin}^{2}\left(\theta_{1}\right) + m_{2} \operatorname{sin}^{2}\left(\theta_{2}\right)\right)}")
        t2ddRHS.scale(0.45)
        t2ddLHS.next_to(qt2, RIGHT)
        t2ddRHS.next_to(t2ddLHS,RIGHT)
        t2dd = VGroup(t2ddLHS, t2ddRHS)
        t2dd.next_to(t1dd, DOWN)
        
        self.play(FadeOutAndShift(qall,UP), FadeOutAndShift(qx, LEFT), FadeOutAndShift(qt1,LEFT), FadeOutAndShift(qt2,LEFT),
                  Transform(xeq,xdd), Transform(t1eq,t1dd), Transform(t2eq,t2dd))
        
        # Linearize that stuff!
        
        xddLRHS = TexMobject(r"\frac{g m_{1} \theta_{1} - g m_{2} \theta_{2} + u}{M}")
        xddLRHS.next_to(xddLHS, RIGHT)
        xddL = VGroup(xddLHS, xddLRHS)
        
        t1ddLRHS = TexMobject(r"\frac{- g m_{2} \theta_{2} + g \left(M + m_{1}\right) \theta_{1} + u}{M l_{1}}")
        t1ddLRHS.next_to(t1ddLHS,RIGHT)
        t1ddL = VGroup(t1ddLHS, t1ddLRHS)
        
        t2ddLRHS = TexMobject(r"\frac{- g m_{1} \theta_{1} + g \left(M + m_{2}\right) \theta_{2} - u}{M l_{2}}")
        t2ddLRHS.next_to(t2ddLHS,RIGHT)
        t2ddL = VGroup(t2ddLHS, t2ddLRHS)
#        self.remove(xdd, t1dd, t2dd)
        self.play(TransformFromCopy(xdd,xddL), TransformFromCopy(t1dd,t1ddL), TransformFromCopy(t2dd,t2ddL),
                  ApplyMethod(xdd.remove), ApplyMethod(t1dd.remove), ApplyMethod(t2dd.remove))
        
        self.wait(1)
        
        
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
        sys = VGroup(cart, p1, p2)
        # Add elements to scene
        self.add(cart, p1, p2)
        self.bring_to_back(rod1, rod2)
        rod1.always_continually_update = True
        print(rod1.points[0]+OUT)
        print(p1.points)
        movement = 0.2*DOWN
#        self.play(ApplyMethod(cart.move_to, self.relVec(cart, movement)), ApplyMethod(p1.move_to, self.relVec(p1, movement)), ApplyMethod(p2.move_to, self.relVec(p2, movement)))
#        self.play(self.moveSys(sys, 0.2*DOWN), Rotate(p1, angle=PI/8, about_point=rod1.points[0]+OUT))
        self.play(self.moveSys(sys, 0.2*DOWN), Rotate(p1, angle=PI/8, about_point=rod1.points[0]+OUT))
#        self.play(FadeIn(m1), FadeIn(m2), FadeIn(rod1), FadeIn(rod2))
#        self.bring_to_back(rod1, rod2)
        
#        self.play(Rotate(p1, angle=PI/8, about_point=rod1.points[0]+OUT))
        self.wait(2)
        
    def getPos(self, obj):
        objpos = [obj.get_x(), obj.get_y(), obj.get_z()]
        return objpos
            
    def relVec(self, obj, pos):
        # generate the vector to move an object to a new position
        return pos+self.getPos(obj)
    
    def rotatePend(self, obj, angle):
        pass
    
    def moveSys(self, sys, x):
        
        return ApplyMethod(sys.move_to, self.relVec(sys, x))
    
#    class Pendulum(Circle, Line, radius, length, mass, angle, position):
#        CONFIG = {
#                "radius": radius,
#                "fill_color": GREY_1,
#                "stroke_color": NEON_GREEN,
#                "fill_opacity": 0.5,
#                }
#        def __init__(self, **kwargs):
#            Circle.__init__(**kwargs)