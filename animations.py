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
        self.energy2L()
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
        self.wait(1)
        
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
        self.wait(1)
        self.play(TransformFromCopy(qall, qx), TransformFromCopy(qall, qt1), TransformFromCopy(qall, qt2))
        self.wait(1)
        
        self.play(Write(xeq), Write(t1eq), Write(t2eq))
        self.wait(1)
        
        
        
        # Equations of motion
#        xdd = TexMobject(r"\ddot{x} = \frac{0.5 g m_{1} \operatorname{sin}\left(2\theta_{1}\right) - 0.5 g m_{2} \operatorname{sin}\left(2\theta_{2}\right) - l_{1} m_{1} \operatorname{sin}\left(\theta_{1}\right) \dot{\theta}_{1}^{2} + l_{2} m_{2} \operatorname{sin}\left(\theta_{2}\right) \dot{\theta}_{2}^{2} + u}{M + m_{1} \operatorname{sin}^{2}\left(\theta_{1}\right) + m_{2} \operatorname{sin}^{2}\left(\theta_{2}\right)}")
#        xdd.scale(0.45)
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
                  ReplacementTransform(xeq,xdd), ReplacementTransform(t1eq,t1dd), ReplacementTransform(t2eq,t2dd))
        self.wait(1)
        
        # Linearize that stuff!
        
        xddL = TexMobject(r"\ddot{x} = \frac{g m_{1} \theta_{1} - g m_{2} \theta_{2} + u}{M}")
        
        t1ddL  = TexMobject(r"\ddot{\theta_1} = \frac{- g m_{2} \theta_{2} + g \left(M + m_{1}\right) \theta_{1} + u}{M l_{1}}")
        
        t2ddL = TexMobject(r"\ddot{\theta_2} = \frac{- g m_{1} \theta_{1} + g \left(M + m_{2}\right) \theta_{2} - u}{M l_{2}}")
        
        t1ddL.center
        t1ddL.shift(1*DOWN)
        xddL.next_to(t1ddL, UP)
        t2ddL.next_to(t1ddL, DOWN)
        
        self.play(ReplacementTransform(xdd,xddL), ReplacementTransform(t1dd,t1ddL), ReplacementTransform(t2dd,t2ddL))
        
        # Break point for the controls section
        self.play(FadeOutAndShift(Lgroup, UP), FadeOutAndShift(Leq3, UP), ApplyMethod(xddL.to_edge, UP), MaintainPositionRelativeTo(t1ddL, xddL), MaintainPositionRelativeTo(t2ddL, xddL))
        
        self.wait(1)

class Controls(Scene):
    def construct(self):
        self.stateSpace()
        self.controllability()
        self.observability()
        self.LQR()
        
    def stateSpace(self):
        xddL = TexMobject(r"\ddot{x} \,", r"= \frac{g m_{1} \theta_{1} - g m_{2} \theta_{2} + u}{M}")
        t1ddL  = TexMobject(r"\ddot{\theta_1} \,", r" = \frac{- g m_{2} \theta_{2} + g \left(M + m_{1}\right) \theta_{1} + u}{M l_{1}}")
        t2ddL = TexMobject(r"\ddot{\theta_2} \,", r" = \frac{- g m_{1} \theta_{1} + g \left(M + m_{2}\right) \theta_{2} - u}{M l_{2}}")
        
        
        xddL.to_edge(UP)
        t1ddL.next_to(xddL, DOWN)
        t2ddL.next_to(t1ddL, DOWN)
        
        for var in [xddL, t1ddL, t2ddL]:
            var[0].scale(1.75)
        
        self.add(xddL, t1ddL, t2ddL)
#        self.play(
#                 xddL[0].scale, 2,
#                 run_time=1,)
        self.play(
                xddL.scale, 0.55,
                xddL.to_corner, UP+LEFT,
                t1ddL.scale, 0.55,
                t1ddL.to_edge, UP,
                t2ddL.scale, 0.55,
                t2ddL.to_corner, UP+RIGHT,
                run_time=1,)
        
        self.wait(1)
        # State vector
        z = TexMobject(r"z = ", r"\left[\begin{matrix} x \\ \dot{x} \\ \theta_1 \\ \dot{\theta_1} \\ \theta_2 \\ \dot{\theta_2} \end{matrix}\right]")
        zd = TexMobject(r"\dot{z} =",  r"\left[\begin{matrix} \dot{x} \\ \ddot{x} \\ \dot{\theta_1} \\ \ddot{\theta_1} \\ \dot{\theta_2} \\ \ddot{\theta_2} \end{matrix}\right]")
        zvec = z[1].copy()
        
        zdexpr = TexMobject(r"\dot{z}", r" =", r"A", r"z", r"+", r"B", r'u')
        deriv = TexMobject(r"\frac{\partial}{\partial t}")
        
        deriv.next_to(z, LEFT)
        self.play(Write(z))
        self.wait(2)
        self.play(Write(deriv))
        self.wait(2)
        self.play(
                ReplacementTransform(z, zd),
                deriv.fade,1)
        self.wait(2)
        self.play(
                zd[1].to_edge, LEFT,
                zd[0].fade, 1,
                run_time=1)
        self.wait(2)
        self.play(
                Write(zdexpr),
                run_time=1,)
        self.wait(2)
        self.play(
                zdexpr.next_to, t1ddL, DOWN,
                run_time=1)
        self.wait(2)
        # Write out A and B
        A = TexMobject(r"= ", r"\left[\begin{matrix}0 & 1 & 0 & 0 & 0 & 0\\0 & 0 & \frac{g m_{1}}{M} & 0 & - \frac{g m_{2}}{M} & 0\\0 & 0 & 0 & 1 & 0 & 0\\0 & 0 & \frac{g \left(M + m_{1}\right)}{M l_{1}} & 0 & - \frac{g m_{2}}{M l_{1}} & 0\\0 & 0 & 0 & 0 & 0 & 1\\0 & 0 & - \frac{g m_{1}}{M l_{2}} & 0 & \frac{g \left(M + m_{2}\right)}{M l_{2}} & 0\end{matrix}\right]")
        
        A.next_to(zd[1], RIGHT)
        zvec.next_to(A, RIGHT)
        self.play(
                TransformFromCopy(zdexpr[2], A),
                TransformFromCopy(zdexpr[3], zvec),
                run_time=1)
        self.wait(2)
        B = TexMobject(r"B = ", r"\left[\begin{matrix}0\\\frac{1}{M}\\0\\\frac{1}{M l_{1}}\\0\\- \frac{1}{M l_{2}}\end{matrix}\right]")
        plus = TexMobject(r"+")
        u = TexMobject(r"u")
        plus.next_to(zvec, RIGHT)
        B[1].next_to(plus, RIGHT)
        u.next_to(B, RIGHT)
        # Full zdot = Az+Bu expression
        self.play(
                TransformFromCopy(zdexpr[4], plus),
                TransformFromCopy(zdexpr[5], B[1]),
                TransformFromCopy(zdexpr[6], u),
                run_time=1)
        
        self.wait(2)
        # Clear top of screen
        zdfullexpr = VGroup(zd[1], A, zvec, plus, B[1], u)
        self.play(
                FadeOutAndShift(xddL, UP),
                FadeOutAndShift(t1ddL, UP),
                FadeOutAndShift(t2ddL, UP),
                FadeOutAndShift(zdexpr, UP),
                zdfullexpr.to_edge, UP,
                run_time=1)
        
        self.wait(2)

        
        # Now for y = C z
        yexpr = TexMobject(r"y", r"=", r"C", r"z")
        yexpr.to_edge(DOWN)
        yexpr.shift(0.75*UP)
#        yexpr.align_to(zd[1], LEFT)
        self.play(Write(yexpr))
        self.wait(2)
        # Reduce the zd = Az+Bu expresion to drive home what z is
        zdsym = TexMobject(r'\dot{z}')
        zsym = TexMobject(r'z')
        
        zdsym.move_to(zd[1])
        zsym.move_to(zvec)
        self.play(
                ReplacementTransform(zd[1],zdsym),
                ReplacementTransform(zvec, zsym),
                plus.next_to, zsym, RIGHT,
                MaintainPositionRelativeTo(B[1], plus),
                MaintainPositionRelativeTo(u, plus))
        # Collapse equation
        self.play(
                zdsym.next_to, A, LEFT,
                zsym.next_to, A, RIGHT,
                MaintainPositionRelativeTo(plus, zsym),
                MaintainPositionRelativeTo(B[1], zsym),
                MaintainPositionRelativeTo(u, zsym))
        self.wait(2)
        yfull = TexMobject(r"y", r"=", r"\left[\begin{matrix}1 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 1 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 1 & 0\end{matrix}\right]", r'z')
        yfull.to_edge(DOWN)
        self.play(
                ReplacementTransform(yexpr, yfull))
        self.wait(2)
        
        # Set up controllablity section
        self.play(
                FadeOutAndShift(zdsym, LEFT),
                FadeOutAndShift(A[0], LEFT),
                zsym.fade, 1,
                FadeOutAndShift(u, RIGHT),
                FadeOutAndShift(yfull, DOWN),
                plus.fade, 1,)
        
        
        self.play(
                A[1].center)
        self.wait(1)
        self.play(B[1].next_to, A[1], RIGHT)
        self.remove(A[1], B[1])
        
    def controllability(self):
        A = TexMobject(r'\left[\begin{matrix}0 & 1 & 0 & 0 & 0 & 0\\0 & 0 & \frac{g m_{1}}{M} & 0 & - \frac{g m_{2}}{M} & 0\\0 & 0 & 0 & 1 & 0 & 0\\0 & 0 & \frac{g \left(M + m_{1}\right)}{M l_{1}} & 0 & - \frac{g m_{2}}{M l_{1}} & 0\\0 & 0 & 0 & 0 & 0 & 1\\0 & 0 & - \frac{g m_{1}}{M l_{2}} & 0 & \frac{g \left(M + m_{2}\right)}{M l_{2}} & 0\end{matrix}\right]')
        B = TexMobject(r'\left[\begin{matrix}0\\\frac{1}{M}\\0\\\frac{1}{M l_{1}}\\0\\- \frac{1}{M l_{2}}\end{matrix}\right]')
        A.center
        B.next_to(A, RIGHT)
        self.add(A, B)
        self.play(A.to_edge, UP,
                  B.to_edge, UP,
                  )
        # Controllability matrix
        Qcexpr = TexMobject(r'Q_c', r'=', r'\left[ B \quad A B \quad A^2 B \quad A^3 B \quad A^4 B \quad A^5 B\right]')
        Qcexpr.next_to(A, DOWN)
        self.play(Write(Qcexpr))
        self.wait(1)
        
        self.play(
                FadeOutAndShift(A, UP),
                FadeOutAndShift(B, UP),
                Qcexpr.to_edge, UP)
        self.wait(1)
#        self.play(Qcexpr.stretch, 1.5, 1)
        
        # Actual Qc
        Qc = TexMobject(r'Q_c', r'=', r'\left[\begin{matrix}0 & \frac{1}{M} & 0 & \frac{g \left(l_{1} m_{2} + l_{2} m_{1}\right)}{M^{2} l_{1} l_{2}} & 0 & \frac{g^{2} \left(l_{1} m_{2} \left(l_{1} \left(M + m_{2}\right) + l_{2} m_{1}\right) + l_{2} m_{1} \left(l_{1} m_{2} + l_{2} \left(M + m_{1}\right)\right)\right)}{M^{3} l_{1}^{2} l_{2}^{2}}\\\frac{1}{M} & 0 & \frac{g \left(l_{1} m_{2} + l_{2} m_{1}\right)}{M^{2} l_{1} l_{2}} & 0 & \frac{g^{2} \left(l_{1} m_{2} \left(l_{1} \left(M + m_{2}\right) + l_{2} m_{1}\right) + l_{2} m_{1} \left(l_{1} m_{2} + l_{2} \left(M + m_{1}\right)\right)\right)}{M^{3} l_{1}^{2} l_{2}^{2}} & 0\\0 & \frac{1}{M l_{1}} & 0 & \frac{g \left(l_{1} m_{2} + l_{2} \left(M + m_{1}\right)\right)}{M^{2} l_{1}^{2} l_{2}} & 0 & \frac{g^{2} \left(l_{1} m_{2} \left(l_{1} \left(M + m_{2}\right) + l_{2} \left(M + m_{1}\right)\right) + l_{2} \left(l_{1} m_{1} m_{2} + l_{2} \left(M + m_{1}\right)^{2}\right)\right)}{M^{3} l_{1}^{3} l_{2}^{2}}\\\frac{1}{M l_{1}} & 0 & \frac{g \left(l_{1} m_{2} + l_{2} \left(M + m_{1}\right)\right)}{M^{2} l_{1}^{2} l_{2}} & 0 & \frac{g^{2} \left(l_{1} m_{2} \left(l_{1} \left(M + m_{2}\right) + l_{2} \left(M + m_{1}\right)\right) + l_{2} \left(l_{1} m_{1} m_{2} + l_{2} \left(M + m_{1}\right)^{2}\right)\right)}{M^{3} l_{1}^{3} l_{2}^{2}} & 0\\0 & - \frac{1}{M l_{2}} & 0 & - \frac{g \left(l_{1} \left(M + m_{2}\right) + l_{2} m_{1}\right)}{M^{2} l_{1} l_{2}^{2}} & 0 & - \frac{g^{2} \left(l_{1} \left(l_{1} \left(M + m_{2}\right)^{2} + l_{2} m_{1} m_{2}\right) + l_{2} m_{1} \left(l_{1} \left(M + m_{2}\right) + l_{2} \left(M + m_{1}\right)\right)\right)}{M^{3} l_{1}^{2} l_{2}^{3}}\\- \frac{1}{M l_{2}} & 0 & - \frac{g \left(l_{1} \left(M + m_{2}\right) + l_{2} m_{1}\right)}{M^{2} l_{1} l_{2}^{2}} & 0 & - \frac{g^{2} \left(l_{1} \left(l_{1} \left(M + m_{2}\right)^{2} + l_{2} m_{1} m_{2}\right) + l_{2} m_{1} \left(l_{1} \left(M + m_{2}\right) + l_{2} \left(M + m_{1}\right)\right)\right)}{M^{3} l_{1}^{2} l_{2}^{3}} & 0\end{matrix}\right]')
        Qc.scale(0.35)
        Qc.stretch(1.5,1)
        self.play(Write(Qc))
        
        # Determinant of Qc
        QcBrace = Brace(Qc[2],DOWN)
        detQc = TexMobject(r'\operatorname{det}Q_c', r'=', r'\frac{g^{6} \left(l_{1} - l_{2}\right)^{2}}{M^{6} l_{1}^{6} l_{2}^{6}}')
        detQc.next_to(QcBrace, DOWN)
        self.play(GrowFromCenter(QcBrace), Write(detQc))
        
        # Controllability condition
        noctrl = TexMobject(r'\operatorname{if}\ ', r'l_1', r'=', r'l_2')
        noctrl.next_to(detQc, LEFT)
        detQc0 = TexMobject(r'0')
        detQc0.next_to(detQc[1],RIGHT)
        self.play(Write(noctrl))
        
        self.play(Transform(detQc[2], detQc0))
        self.remove(detQc[2])
        noctrlgroup = VGroup(noctrl, detQc0, detQc[0], detQc[1])
        
        self.play(
                FadeOutAndShift(Qc, UP),
                FadeOutAndShift(QcBrace, UP),
                FadeOutAndShift(Qcexpr, UP),
                noctrlgroup.center)
        
        # State explicitly that it's uncontrollable
#        ctrlcopy = noctrlgroup.copy()
        noctrl2 = TexMobject(r'\operatorname{if}\ ', r'\operatorname{det}Q_c', r'=', r'0')
        self.play(noctrlgroup.shift,1*UP)
        noctrl2.next_to(noctrlgroup,DOWN)
        self.play(TransformFromCopy(noctrlgroup, noctrl2))
        ctrlstatement = TextMobject(r'The system is uncontrollable')
        ctrlstatement.next_to(noctrl2, DOWN)
        self.play(Write(ctrlstatement))
        
        self.play(FadeOut(ctrlstatement),
                  FadeOut(noctrl2),
                  FadeOut(noctrlgroup))
        
        
    def observability(self):
        ##
        # Observability
        ##
        AandC = TexMobject(r'A = ', r'\left[\begin{matrix}0 & 1 & 0 & 0 & 0 & 0\\0 & 0 & \frac{g m_{1}}{M} & 0 & - \frac{g m_{2}}{M} & 0\\0 & 0 & 0 & 1 & 0 & 0\\0 & 0 & \frac{g \left(M + m_{1}\right)}{M l_{1}} & 0 & - \frac{g m_{2}}{M l_{1}} & 0\\0 & 0 & 0 & 0 & 0 & 1\\0 & 0 & - \frac{g m_{1}}{M l_{2}} & 0 & \frac{g \left(M + m_{2}\right)}{M l_{2}} & 0\end{matrix}\right]',
                           r'C = ', r'\left[\begin{matrix}1 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 1 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 1 & 0\end{matrix}\right]')
        AandC.scale(0.75)
        self.play(Write(AandC))
        self.play(AandC.to_edge, UP)
        Qoexpr = TexMobject(r'Q_o', r'=', r'\left[ C \quad C A \quad C A^2 \quad C A^3 \quad C A^4 \quad C A^5 \right]')
        Qoexpr.next_to(AandC, DOWN)
        self.play(Write(Qoexpr))
        self.wait(0.5)
        
        Qo = TexMobject(r'Q_0',r'=', r' \left[\begin{matrix}1 & 0 & 0 & 0 & 0 & 0\\0 & 1 & 0 & 0 & 0 & 0\\0 & 0 & \frac{g m_{1}}{M} & 0 & \frac{g^{2} m_{1} \left(l_{1} m_{2} + l_{2} \left(M + m_{1}\right)\right)}{M^{2} l_{1} l_{2}} & 0\\0 & 0 & 0 & \frac{g m_{1}}{M} & 0 & \frac{g^{2} m_{1} \left(l_{1} m_{2} + l_{2} \left(M + m_{1}\right)\right)}{M^{2} l_{1} l_{2}}\\0 & 0 & - \frac{g m_{2}}{M} & 0 & - \frac{g^{2} m_{2} \left(l_{1} \left(M + m_{2}\right) + l_{2} m_{1}\right)}{M^{2} l_{1} l_{2}} & 0\\0 & 0 & 0 & - \frac{g m_{2}}{M} & 0 & - \frac{g^{2} m_{2} \left(l_{1} \left(M + m_{2}\right) + l_{2} m_{1}\right)}{M^{2} l_{1} l_{2}}\end{matrix}\right]')
        Qo.scale(0.75)
        self.play(
                FadeOutAndShift(AandC, UP),
                Qoexpr.to_edge, UP,
                TransformFromCopy(Qoexpr, Qo))
        self.wait(1)

        
        # Determinant of Qo
        QoBrace = Brace(Qo[2],DOWN)
        detQo = TexMobject(r'\operatorname{det}Q_o', r'=', r'\frac{g^{6} m_{1}^{2} m_{2}^{2} \left(l_{1} - l_{2}\right)^{2}}{M^{4} l_{1}^{2} l_{2}^{2}}')
        detQo.next_to(QoBrace, DOWN)
        self.play(GrowFromCenter(QoBrace), Write(detQo))
        
        # observeability condition
        noobs = TexMobject(r'\operatorname{if}\ ', r'l_1', r'=', r'l_2')
        noobs.next_to(detQo, LEFT)
        detQo0 = TexMobject(r'0')
        detQo0.next_to(detQo[1],RIGHT)
        self.play(Write(noobs))
        
        self.play(Transform(detQo[2], detQo0))
        self.remove(detQo[2])
        noobsgroup = VGroup(noobs, detQo0, detQo[0], detQo[1])
        
        self.play(
                FadeOutAndShift(Qo, UP),
                FadeOutAndShift(QoBrace, UP),
                FadeOutAndShift(Qoexpr, UP),
                noobsgroup.center)
        
        # State explicitly that it's unobserveable
#        obscopy = noobsgroup.copy()
        noobs2 = TexMobject(r'\operatorname{if}\ ', r'\operatorname{det}Q_o', r'=', r'0')
        self.play(noobsgroup.shift,1*UP)
        noobs2.next_to(noobsgroup,DOWN)
        self.play(TransformFromCopy(noobsgroup, noobs2))
        obsstatement = TextMobject(r'The system is unobserveable')
        obsstatement.next_to(noobs2, DOWN)
        self.play(Write(obsstatement))
        self.play(FadeOut(obsstatement),
                  FadeOut(noobs2),
                  FadeOut(noobsgroup))
        
    def LQR(self):
        lqrexpr = TexMobject(r'J = ', r'\int_0^\infty (  ', r'z^T', r'Q' ,r'z', r'+', r'u^T', r'R', ' u', r' ) dt')
        lqrexpr.to_edge(UP)
        self.play(Write(lqrexpr))
        
        
        z = TexMobject(r"z = ", r"\left[\begin{matrix} x \\ \dot{x} \\ \theta_1 \\ \dot{\theta_1} \\ \theta_2 \\ \dot{\theta_2} \end{matrix}\right]")
        z.to_edge(LEFT)
        self.play(TransformFromCopy(lqrexpr[4], z))
        
        
        Q = TexMobject(r'Q', r'=', r'\rho', r'\left[\begin{matrix}0.1 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 1 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 1 & 0\\0 & 0 & 0 & 0 & 0 & 0\end{matrix}\right]')
        self.play(TransformFromCopy(lqrexpr[3], Q))
        
        R = TexMobject(r'R', '=', r'\mu')
        R.to_edge(RIGHT)
        self.play(TransformFromCopy(lqrexpr[7], R))
        
        rhoval = TexMobject(r'\rho' ,'=' , '100')
        Q2 = TexMobject(r'100', r'\left[\begin{matrix}10 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 100 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 100 & 0\\0 & 0 & 0 & 0 & 0 & 0\end{matrix}\right]')
        Q2[0].move_to(Q[2])
        Q2[0].shift(0.075*UP+0.12*LEFT)
        
        muval = TexMobject(r'\mu ', '=', '1')
        muval.next_to(rhoval, RIGHT)
        muval.shift(0.2*RIGHT)
        valsgroup = VGroup(rhoval, muval)
        valsgroup.next_to(Q, DOWN)
        
        self.play(Write(rhoval))
        self.play(Transform(rhoval,Q2[0]),
                  FadeOut(Q[2]),
                  Q[0].shift, 0.3*LEFT,
                  Q[1].shift, 0.3*LEFT,)
        
        Q2[1].next_to(Q[1])
        
        self.play(Transform(Q2[0], Q2[1]),
                  FadeOut(Q[3]),
                  FadeOut(rhoval))
        self.wait(1)
        
        self.play(Write(muval))
        
        Rval = TexMobject(r'1')
        Rval.next_to(R[1], RIGHT)
        self.play(Transform(muval,Rval),
                  FadeOut(R[2]))
        
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