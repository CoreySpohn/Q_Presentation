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
        Texpr = TexMobject(r'KE', '=', r'\sum \frac{1}{2} m_i \vec{v_i} \cdot \vec{v_i}')
        Texpr.to_edge(UP)
        self.play(Write(Texpr))
        Teq = TexMobject(r'KE', '=', r'\frac{1}{2} M \dot{x}^{2}',
                         r'+ \frac{1}{2} m_{1} \left(l_{1}^{2} \dot{\theta}_{1}^{2} - 2 l_{1} \operatorname{cos}\left(\theta_{1}\right) \dot{\theta}_{1} \dot{x} + \dot{x}^{2}\right)',
                         r'+ \frac{1}{2} m_{2} \left(l_{2}^{2} \dot{\theta}_{2}^{2} + 2 l_{2} \operatorname{cos}\left(\theta_{2}\right) \dot{\theta}_{2} \dot{x} + \dot{x}^{2}\right)')
#        Teq[2].shift(5*RIGHT)
        
        Teq[2].center()
        Teq[2].shift(3*LEFT)
        Teq[1].next_to(Teq[2], LEFT)
        Teq[0].next_to(Teq[1], LEFT)
        Teq[3].next_to(Teq[2], DOWN)
        Teq[3].align_to(Teq[2], LEFT)
        Teq[4].next_to(Teq[3], DOWN)
        Teq[4].align_to(Teq[3], LEFT)
        
        self.play(Write(Teq[0]),
                  Write(Teq[1]),
                  Write(Teq[2]))
        self.wait(0.5)
        self.play(Write(Teq[3]))
        self.wait(0.5)
        self.play(Write(Teq[4]))
        self.wait(0.5)
        self.play(FadeOutAndShift(Texpr, UP),
                  Teq.to_edge, UP,
                  Teq.shift, 1*DOWN)
        # Potential Energy
        Uexpr = TexMobject(r'PE', '=', r'\sum g m_i h_i')
        Uexpr.next_to(Teq, DOWN)
        self.play(Write(Uexpr))
        self.wait(0.5)
        Ueq = TexMobject(r"PE =", r" g \left(l_{1} m_{1} \operatorname{cos}\left(\theta_{1}\right) + l_{2} m_{2} \operatorname{cos}\left(\theta_{2}\right)\right)")
        Ueq.move_to(Uexpr)
        Ueq.align_to(Teq, LEFT)
        self.play(ReplacementTransform(Uexpr, Ueq))
        self.wait(0.5)
        
        # Lagrangian symbolic
        Lexpr = TexMobject(r"\mathcal{L} =", r'KE - PE')
        Lexpr.to_edge(UP)
        self.play(Write(Lexpr))
        self.wait(0.5)
        # Lagrangian full
        Leq = TexMobject(r"\mathcal{L} =",
                         r"\frac{1}{2} M \dot{x}^{2} + \frac{1}{2} m_{1} \left(l_{1}^{2} \dot{\theta}_{1}^{2} - 2 l_{1} \operatorname{cos}\left(\theta_{1}\right) \dot{\theta}_{1} \dot{x} + \dot{x}^{2}\right)",
                         r"+ \frac{1}{2} m_{2} \left(l_{2}^{2} \dot{\theta}_{2}^{2} + 2 l_{2} \operatorname{cos}\left(\theta_{2}\right) \dot{\theta}_{2} \dot{x} + \dot{x}^{2}\right)",
                         r"-g \left(l_{1} m_{1} \operatorname{cos}\left(\theta_{1}\right) + l_{2} m_{2} \operatorname{cos}\left(\theta_{2}\right)\right)")
        
        Leq[0].move_to(Lexpr)
        Leq[0].align_to(Teq[0], LEFT)
        Leq[1].next_to(Leq[0], RIGHT)
        Leq[2].next_to(Leq[1], DOWN)
        Leq[2].align_to(Leq[1], LEFT)
        Leq[3].next_to(Leq[2], DOWN)
        Leq[3].align_to(Leq[2], LEFT)
        Leq.center()
        Leq.to_edge(UP)
        
        self.play(ReplacementTransform(Lexpr, Leq[0]))
        self.wait(0.5)
        self.play(ReplacementTransform(Teq[3], Leq[1]),
                  ReplacementTransform(Teq[4], Leq[2]),
                  FadeOut(Teq[2]),
                  FadeOut(Teq[0]),
                  FadeOut(Teq[1]))
        self.wait(0.5)
        self.play(ReplacementTransform(Ueq, Leq[3]))
        
#        self.remove(Ueq)
        self.wait(0.5)
        self.remove(Leq, Teq, Texpr, Ueq, Lexpr, Uexpr, Leq[0], Leq[3], Leq[1], Leq[2])
        
    def l2EOM(self):
        # Lagrange equation 
        Leq = TexMobject(r"\mathcal{L} =",
                         r"\frac{1}{2} M \dot{x}^{2} + \frac{1}{2} m_{1} \left(l_{1}^{2} \dot{\theta}_{1}^{2} - 2 l_{1} \operatorname{cos}\left(\theta_{1}\right) \dot{\theta}_{1} \dot{x} + \dot{x}^{2}\right)",
                         r"+ \frac{1}{2} m_{2} \left(l_{2}^{2} \dot{\theta}_{2}^{2} + 2 l_{2} \operatorname{cos}\left(\theta_{2}\right) \dot{\theta}_{2} \dot{x} + \dot{x}^{2}\right)",
                         r"-g \left(l_{1} m_{1} \operatorname{cos}\left(\theta_{1}\right) + l_{2} m_{2} \operatorname{cos}\left(\theta_{2}\right)\right)")
        

        Leq[1].next_to(Leq[0], RIGHT)
        Leq[2].next_to(Leq[1], DOWN)
        Leq[2].align_to(Leq[1], LEFT)
        Leq[3].next_to(Leq[2], DOWN)
        Leq[3].align_to(Leq[2], LEFT)
        Leq.center()
        Leq.to_edge(UP)
        self.add(Leq)
        # The Lagrange equations of motion
        Leq3 = TexMobject(r" \frac{d}{dt}\left(\frac{\partial \mathcal{L}}{\partial \dot{q}}\right) -\frac{\partial \mathcal{L}}{\partial q} = \sum \frac{\partial \vec{r}}{\partial q} \cdot \vec{F}")
        Leq3.next_to(Leq, DOWN)
        self.play(Write(Leq3))
        # Scale and move
        self.play(Leq.scale,0.75,
                  Leq.to_corner,UP+LEFT,
                  Leq3.scale,0.75,
                  Leq3.to_corner,UP+RIGHT)
        self.wait(1)
        
        # Generalized coordinates
        qall = TexMobject(r"q=(x, \theta_1, \theta_2)")
        qx = TexMobject(r"q=x :")
        qt1 = TexMobject(r"q=\theta_1 :")
        qt2 = TexMobject(r"q=\theta_2 :")
        
        qall.next_to(Leq,DOWN)
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
        self.wait(0.5)
        self.play(TransformFromCopy(qall, qx))
        self.wait(0.5)
        self.play(Write(xeq))
        self.wait(0.5)
        self.play(TransformFromCopy(qall, qt1))
        self.wait(0.5)
        self.play(Write(t1eq))
        self.wait(0.5)
        self.play(TransformFromCopy(qall, qt2))
        self.wait(0.5)
        self.play(Write(t2eq))
        self.wait(0.5)
        
        
        
        # Equations of motion

        xdd = TexMobject(r"\ddot{x} = ", 
                         r"\frac{0.5 g m_{1} \operatorname{sin}\left(2\theta_{1}\right) - 0.5 g m_{2} \operatorname{sin}\left(2\theta_{2}\right) - l_{1} m_{1} \operatorname{sin}\left(\theta_{1}\right) \dot{\theta}_{1}^{2} + l_{2} m_{2} \operatorname{sin}\left(\theta_{2}\right) \dot{\theta}_{2}^{2} + u}{M + m_{1} \operatorname{sin}^{2}\left(\theta_{1}\right) + m_{2} \operatorname{sin}^{2}\left(\theta_{2}\right)}")

        t1dd = TexMobject(r"\ddot{\theta_1} =",
                         r"\frac{- 0.25 g m_{2} \left(- \operatorname{sin}\left(\theta_{1} - 2 \theta_{2}\right) + \operatorname{sin}\left(\theta_{1} + 2 \theta_{2}\right)\right) + g \left(M + m_{1} + m_{2} \operatorname{sin}^{2}\left(\theta_{2}\right)\right) \operatorname{sin}\left(\theta_{1}\right) + \left(- l_{1} m_{1} \operatorname{sin}\left(\theta_{1}\right) \dot{\theta}_{1}^{2} + l_{2} m_{2} \operatorname{sin}\left(\theta_{2}\right) \dot{\theta}_{2}^{2} + u\right) \operatorname{cos}\left(\theta_{1}\right)}{l_{1} \left(M + m_{1} \operatorname{sin}^{2}\left(\theta_{1}\right) + m_{2} \operatorname{sin}^{2}\left(\theta_{2}\right)\right)}")

        t2dd = TexMobject(r"\ddot{\theta_2} =",
                         r"\frac{- 0.25 g m_{1} \left(\operatorname{sin}\left(2 \theta_{1} - \theta_{2}\right) + \operatorname{sin}\left(2 \theta_{1} + \theta_{2}\right)\right) + g \left(M + m_{1} \operatorname{sin}^{2}\left(\theta_{1}\right) + m_{2}\right) \operatorname{sin}\left(\theta_{2}\right) - \left(- l_{1} m_{1} \operatorname{sin}\left(\theta_{1}\right) \dot{\theta}_{1}^{2} + l_{2} m_{2} \operatorname{sin}\left(\theta_{2}\right) \dot{\theta}_{2}^{2} + u\right) \operatorname{cos}\left(\theta_{2}\right)}{l_{2} \left(M + m_{1} \operatorname{sin}^{2}\left(\theta_{1}\right) + m_{2} \operatorname{sin}^{2}\left(\theta_{2}\right)\right)}")
        
        for var in [xdd, t1dd, t2dd]:
            var[0].scale(1)
            var[1].scale(0.4)
        
        xdd[0].next_to(qx, RIGHT)
        xdd[0].to_edge(LEFT)
        xdd[1].next_to(xdd[0],RIGHT)
        t1dd[0].next_to(xdd[0], DOWN)
        t1dd[1].next_to(t1dd[0],RIGHT)
        t1dd.shift(0.75*DOWN)
        t2dd[0].next_to(t1dd[0], DOWN)
        t2dd[1].next_to(t2dd[0],RIGHT)
        t2dd.shift(0.75*DOWN)
        
        self.play(FadeOutAndShift(qall,UP), FadeOutAndShift(qx, LEFT), FadeOutAndShift(qt1,LEFT), FadeOutAndShift(qt2,LEFT))
        self.play(ReplacementTransform(xeq,xdd))
        self.wait(0.5)
        self.play(ReplacementTransform(t1eq,t1dd))
        self.wait(0.5)
        self.play(ReplacementTransform(t2eq,t2dd))
        self.wait(0.5)
        
        # Linearize that stuff!
        
        xddL = TexMobject(r"\ddot{x} =", r"\frac{g m_{1} \theta_{1} - g m_{2} \theta_{2} + u}{M}")
        t1ddL  = TexMobject(r"\ddot{\theta_1} =", r"\frac{- g m_{2} \theta_{2} + g \left(M + m_{1}\right) \theta_{1} + u}{M l_{1}}")
        t2ddL = TexMobject(r"\ddot{\theta_2} =", r"\frac{- g m_{1} \theta_{1} + g \left(M + m_{2}\right) \theta_{2} - u}{M l_{2}}")
        
        t1ddL.center
        t1ddL.shift(1*DOWN)
        xddL.next_to(t1ddL, UP)
        t2ddL.next_to(t1ddL, DOWN)
        
        self.play(ReplacementTransform(xdd,xddL))
        self.play(ReplacementTransform(t1dd,t1ddL))
        self.play(ReplacementTransform(t2dd,t2ddL))
        self.wait(0.5)
        # Break point for the controls section
        self.play(FadeOutAndShift(Leq, UP), FadeOutAndShift(Leq3, UP), 
                  ApplyMethod(xddL.to_edge, UP), MaintainPositionRelativeTo(t1ddL, xddL),
                  MaintainPositionRelativeTo(t2ddL, xddL))
        
        self.wait(2)

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
        
        self.add(xddL, t1ddL, t2ddL)

        self.play(
                xddL.scale, 0.55,
                xddL.to_corner, UP+LEFT,
                t1ddL.scale, 0.55,
                t1ddL.to_edge, UP,
                t2ddL.scale, 0.55,
                t2ddL.to_corner, UP+RIGHT,
                run_time=1,)
        
        self.wait(0.5)
        # State vector
        z = TexMobject(r"z = ", r"\left[\begin{matrix} x \\ \dot{x} \\ \theta_1 \\ \dot{\theta_1} \\ \theta_2 \\ \dot{\theta_2} \end{matrix}\right]")
        zd = TexMobject(r"\dot{z} =",  r"\left[\begin{matrix} \dot{x} \\ \ddot{x} \\ \dot{\theta_1} \\ \ddot{\theta_1} \\ \dot{\theta_2} \\ \ddot{\theta_2} \end{matrix}\right]")
        zvec = z[1].copy()
        
        zdexpr = TexMobject(r"\dot{z}", r" =", r"A", r"z", r"+", r"B", r'u')
        deriv = TexMobject(r"\frac{\partial}{\partial t}")
        
        deriv.next_to(z, LEFT)
        self.play(Write(z))
        self.wait(0.5)
        self.play(Write(deriv))
        self.wait(0.5)
        self.play(
                ReplacementTransform(z, zd),
                deriv.fade,1)
        self.wait(0.5)
        self.play(
                zd[1].to_edge, LEFT,
                zd[0].fade, 1,
                run_time=1)
        self.wait(0.5)
        self.play(
                Write(zdexpr),
                run_time=1,)
        self.wait(2)
        self.play(
                zdexpr.next_to, t1ddL, DOWN,
                run_time=1)
        self.wait(0.5)
        # Write out A and B
        A = TexMobject(r"= ", r"\left[\begin{matrix}0 & 1 & 0 & 0 & 0 & 0\\0 & 0 & \frac{g m_{1}}{M} & 0 & - \frac{g m_{2}}{M} & 0\\0 & 0 & 0 & 1 & 0 & 0\\0 & 0 & \frac{g \left(M + m_{1}\right)}{M l_{1}} & 0 & - \frac{g m_{2}}{M l_{1}} & 0\\0 & 0 & 0 & 0 & 0 & 1\\0 & 0 & - \frac{g m_{1}}{M l_{2}} & 0 & \frac{g \left(M + m_{2}\right)}{M l_{2}} & 0\end{matrix}\right]")
        
        A.next_to(zd[1], RIGHT)
        zvec.next_to(A, RIGHT)
        self.play(
                TransformFromCopy(zdexpr[2], A),
                TransformFromCopy(zdexpr[3], zvec),
                run_time=1)
        self.wait(0.5)
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
        
        self.wait(1)
        # Clear top of screen
        zdfullexpr = VGroup(zd[1], A, zvec, plus, B[1], u)
        self.play(
                FadeOutAndShift(xddL, UP),
                FadeOutAndShift(t1ddL, UP),
                FadeOutAndShift(t2ddL, UP),
                FadeOutAndShift(zdexpr, UP),
                zdfullexpr.to_edge, UP,
                run_time=1)
        
        self.wait(1)

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
        
        # Now for y = C z
        yexpr = TexMobject(r"y", r"=", r"C", r"z")
        yexpr.to_edge(DOWN)
        yexpr.shift(0.75*UP)
#        yexpr.align_to(zd[1], LEFT)
        self.play(Write(yexpr))
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
        
       
        
        A2 = TexMobject(r"A = ", r"\left[\begin{matrix}0 & 1 & 0 & 0 & 0 & 0\\0 & 0 & \frac{g m_{1}}{M} & 0 & - \frac{g m_{2}}{M} & 0\\0 & 0 & 0 & 1 & 0 & 0\\0 & 0 & \frac{g \left(M + m_{1}\right)}{M l_{1}} & 0 & - \frac{g m_{2}}{M l_{1}} & 0\\0 & 0 & 0 & 0 & 0 & 1\\0 & 0 & - \frac{g m_{1}}{M l_{2}} & 0 & \frac{g \left(M + m_{2}\right)}{M l_{2}} & 0\end{matrix}\right]")
        B2 = TexMobject(r"B = ", r"\left[\begin{matrix}0\\\frac{1}{M}\\0\\\frac{1}{M l_{1}}\\0\\- \frac{1}{M l_{2}}\end{matrix}\right]")
        B2.next_to(A2, RIGHT)
        newmatgroup = VGroup(A2, B2)
        newmatgroup.center()
        newmatgroup.to_edge(UP)
        self.play(Transform(A[1], A2),
                  Transform(B[1], B2))
        self.remove(A2, B2, A[1], B[1])
        
    def controllability(self):
        A = TexMobject(r'A = ', r'\left[\begin{matrix}0 & 1 & 0 & 0 & 0 & 0\\0 & 0 & \frac{g m_{1}}{M} & 0 & - \frac{g m_{2}}{M} & 0\\0 & 0 & 0 & 1 & 0 & 0\\0 & 0 & \frac{g \left(M + m_{1}\right)}{M l_{1}} & 0 & - \frac{g m_{2}}{M l_{1}} & 0\\0 & 0 & 0 & 0 & 0 & 1\\0 & 0 & - \frac{g m_{1}}{M l_{2}} & 0 & \frac{g \left(M + m_{2}\right)}{M l_{2}} & 0\end{matrix}\right]')
        B = TexMobject(r'B = ', r'\left[\begin{matrix}0\\\frac{1}{M}\\0\\\frac{1}{M l_{1}}\\0\\- \frac{1}{M l_{2}}\end{matrix}\right]')
        B.next_to(A, RIGHT)
        newmatgroup = VGroup(A, B)
        newmatgroup.center()
        newmatgroup.to_edge(UP)
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
        self.wait(0.5)
        # Determinant of Qc
        QcBrace = Brace(Qc[2],DOWN)
        detQc = TexMobject(r'\operatorname{det}Q_c', r'=', r'\frac{g^{6} \left(l_{1} - l_{2}\right)^{2}}{M^{6} l_{1}^{6} l_{2}^{6}}')
        detQc.next_to(QcBrace, DOWN)
        self.play(GrowFromCenter(QcBrace), Write(detQc))
        self.wait(0.5)
        
        # Controllability condition
        noctrl = TexMobject(r'\operatorname{if}\ ', r'l_1', r'=', r'l_2')
        noctrl.next_to(detQc, LEFT)
        detQc0 = TexMobject(r'0')
        detQc0.next_to(detQc[1],RIGHT)
        self.play(Write(noctrl))
        self.wait(0.5)
        
        self.play(Transform(detQc[2], detQc0))
        self.wait(0.5)
        self.remove(detQc[2])
        noctrlgroup = VGroup(noctrl, detQc0, detQc[0], detQc[1])
        
        self.play(
                FadeOutAndShift(Qc, UP),
                FadeOutAndShift(QcBrace, UP),
                FadeOutAndShift(Qcexpr, UP),
                noctrlgroup.center)
        self.wait(0.5)
        # State explicitly that it's uncontrollable
#        ctrlcopy = noctrlgroup.copy()
        noctrl2 = TexMobject(r'\operatorname{if}\ ', r'\operatorname{det}Q_c', r'=', r'0')
#        self.play(noctrlgroup.shift,1*UP)
        noctrl2.next_to(noctrlgroup,DOWN)
        self.play(TransformFromCopy(noctrlgroup, noctrl2))
        ctrlstatement = TextMobject(r'The system is uncontrollable')
        ctrlstatement.next_to(noctrl2, DOWN)
        self.play(Write(ctrlstatement))
        self.wait(0.5)
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
        AandC.to_edge(UP)
        self.play(Write(AandC))
#        self.play(AandC.to_edge, UP)
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
        self.wait(0.5)
        # observeability condition
        noobs = TexMobject(r'\operatorname{if}\ ', r'l_1', r'=', r'l_2')
        noobs.next_to(detQo, LEFT)
        detQo0 = TexMobject(r'0')
        detQo0.next_to(detQo[1],RIGHT)
        self.play(Write(noobs))
        self.wait(0.5)
        self.play(Transform(detQo[2], detQo0))
        self.wait(0.5)
        self.remove(detQo[2])
        noobsgroup = VGroup(noobs, detQo0, detQo[0], detQo[1])
        
        self.play(
                FadeOutAndShift(Qo, UP),
                FadeOutAndShift(QoBrace, UP),
                FadeOutAndShift(Qoexpr, UP),
                noobsgroup.center)
        self.wait(0.5)
        # State explicitly that it's unobserveable
#        obscopy = noobsgroup.copy()
        noobs2 = TexMobject(r'\operatorname{if}\ ', r'\operatorname{det}Q_o', r'=', r'0')
#        self.play(noobsgroup.shift,1*UP)
        noobs2.next_to(noobsgroup,DOWN)
        self.play(TransformFromCopy(noobsgroup, noobs2))
        self.wait(0.5)
        obsstatement = TextMobject(r'The system is unobservable')
        obsstatement.next_to(noobs2, DOWN)
        self.play(Write(obsstatement))
        self.wait(0.5)
        self.play(FadeOut(obsstatement),
                  FadeOut(noobs2),
                  FadeOut(noobsgroup))
        self.wait(0.5)
    def LQR(self):
        # State feedback
        nofeedback = TexMobject(r'\dot{z}=', r'Az',r'+', 'B u')
        nofeedback.to_edge(UP)
        self.play(Write(nofeedback))
        self.wait(0.5)
        uexpr = TexMobject(r'u', '=', '-', 'K z')
        uexpr.next_to(nofeedback, DOWN)
        self.play(Write(uexpr))
        
        feedback = TexMobject(r'\dot{z} = (A-B K)z')
        feedback.move_to(nofeedback)
        
        self.play(ReplacementTransform(nofeedback, feedback))
        self.wait(0.5)
        self.play(feedback.to_corner, UP+LEFT,
                  uexpr.to_corner, UP+RIGHT)
        lqrexpr = TexMobject(r'J = ', r'\int_0^\infty (  ', r'z^T', r'Q' ,r'z', r'+', r'u^T', r'R', ' u', r' ) dt')
        lqrexpr.to_edge(UP)
        self.play(Write(lqrexpr))
        self.wait(0.5)
        
        z = TexMobject(r"z = ", r"\left[\begin{matrix} x \\ \dot{x} \\ \theta_1 \\ \dot{\theta_1} \\ \theta_2 \\ \dot{\theta_2} \end{matrix}\right]")
        z.to_edge(LEFT)
        self.play(TransformFromCopy(lqrexpr[4], z))
        self.wait(0.5)
        
        Q = TexMobject(r'Q', r'=', r'\rho', r'\left[\begin{matrix}0.01 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 1 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 1 & 0\\0 & 0 & 0 & 0 & 0 & 0\end{matrix}\right]')
        self.play(TransformFromCopy(lqrexpr[3], Q))
        self.wait(0.5)
        R = TexMobject(r'R', '=', r'\mu')
        R.to_edge(RIGHT)
        self.play(TransformFromCopy(lqrexpr[7], R))
        
        rhoval = TexMobject(r'\rho' ,'=' , '100')
        Q2 = TexMobject(r'100', r'\left[\begin{matrix}1 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 100 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 0 & 0\\0 & 0 & 0 & 0 & 100 & 0\\0 & 0 & 0 & 0 & 0 & 0\end{matrix}\right]')
        Q2[0].move_to(Q[2])
        Q2[0].shift(0.075*UP+0.12*LEFT)
        
        muval = TexMobject(r'\mu ', '=', '1')
        muval.next_to(rhoval, RIGHT)
        muval.shift(0.2*RIGHT)
        valsgroup = VGroup(rhoval, muval)
        valsgroup.next_to(Q, DOWN)
        self.wait(0.5)
        self.play(Write(rhoval))
        self.wait(0.5)
        self.play(Transform(rhoval,Q2[0]),
                  FadeOut(Q[2]),
                  Q[0].shift, 0.3*LEFT,
                  Q[1].shift, 0.3*LEFT,)
        
        Q2[1].next_to(Q[1])
        self.wait(0.5)
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
        
        self.play(FadeOut(Q2[0]), FadeOut(Q[1]), FadeOut(Q[0]), FadeOut(R[0]), FadeOut(R[1]), FadeOut(muval), FadeOut(z))
        
        Pexpr = TexMobject(r'0 = P A + A^T P - P B R^{-1}B^T P + Q')
        Pexpr.center()
        ueq = TexMobject(r'u=-', r'R^{-1} B^T P', 'z')
        ueq.next_to(Pexpr, DOWN)
        Kbrace = Brace(ueq[1])
        K = TexMobject(r'K')
        K.next_to(Kbrace, DOWN)
        self.play(Write(Pexpr))
        self.wait(0.5)
        self.play(Write(ueq))
        self.wait(0.5)
        self.play(GrowFromCenter(Kbrace), Write(K))
        self.wait(0.5)
        
class CartDemo(Scene):
    def construct(self):
#        self.simpleMotion()
        self.initialConditionPic()
        
    def simpleMotion(self):
        m1radius = 0.25
        m2radius = 0.25
        
        # Lengths, should also add angle
        l1 = 4
        l2 = 4
        
        # Angles (rad)
        t1 = PI/7 
        t2 = PI/7 
        
        # Unit vectors
        er1 = -np.sin(t1)*RIGHT+np.cos(t1)*UP
        et1 = -np.cos(t1)*RIGHT-np.sin(t1)*UP
        
        er2 = np.sin(t2)*RIGHT+np.cos(t2)*UP
        et2 = np.cos(t2)*RIGHT-np.sin(t2)*UP
        
        # Cart properties
        r_M_O = [0, -1] # Position of cart wrt the origin
        d_M = [0.75, 5] # Cart height, width
        
        # Wheel properties
        wheelradius = 0.4
        r_w1_O = [-3*wheelradius+r_M_O[0], -0.25*d_M[0]-wheelradius+r_M_O[1]]
        r_w2_O = [3*wheelradius+r_M_O[0], -0.25*d_M[0]-wheelradius+r_M_O[1]]
        
        # Ground properties
        r_G = [0,-0.5*d_M[0]-1.64*wheelradius+r_M_O[1],0]
        
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
        rod1 = Line(np.array(r_C1_O), np.array(r_m1_O), stroke_color=GREY_1)
        rod2 = Line(np.array(r_C2_O), np.array(r_m2_O), stroke_color=GREY_1)
        
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
        
        # text 
        Mtext = TexMobject(r'M')
        Mtext.move_to(M)
        
        m1text_1 = TexMobject(r'm')
        m1text_1.next_to(m1, UP)
        l1text_1 = TexMobject(r'l')
        l1text_1.next_to(rod1, LEFT)
        l1text_1.shift(0.75*RIGHT)
        
        m2text_1 = TexMobject(r'm')
        m2text_1.next_to(m2, UP)
        l2text_1 = TexMobject(r'l')
        l2text_1.next_to(rod2, RIGHT)
        l2text_1.shift(0.75*LEFT)
        
        # Create groups to move together
        cart = VGroup(M, w1, w2)
        
        p1 = VGroup(m1, rod1)
        p2 = VGroup(m2, rod2)
        sys = VGroup(cart, p1, p2)
        
        # Add elements to scene
        # Create ground
        ground = Line(np.array(r_G+6*LEFT), np.array(r_G+6*RIGHT))
        self.play(Write(ground),run_time=0.3)
        self.play(DrawBorderThenFill(cart),Write(Mtext))
        self.wait(1)
        self.play(DrawBorderThenFill(m1), DrawBorderThenFill(rod1),
                  DrawBorderThenFill(m2), DrawBorderThenFill(rod2),run_time=0.7)
        self.wait(1)
        
        self.play(Write(m1text_1), Write(m2text_1))
        self.play(Write(l1text_1), Write(l2text_1))
        
        ##
        # Angles and position x
        ##
        # theta1
        t1refline = DashedLine(np.array(r_C1_O), np.array(r_C1_O+7*UP), stroke_color=WHITE)
        t1arcpoint1 = r_C1_O + 0.5*l1*UP
        t1arcpoint2 = r_C1_O + 0.5*l1*er1
        t1arc = ArcBetweenPoints(t1arcpoint1, t1arcpoint2)
        t1text = TexMobject(r'\theta_1')
        t1text.next_to(t1arc, UP)
        t1text.scale(1.25)
        self.play(Write(t1refline),
                  run_time=0.5)
        self.play(Write(t1arc),
                  run_time=0.5)
        self.play(Write(t1text))
        
        # theta2
        t2refline = DashedLine(np.array(r_C2_O), np.array(r_C2_O+7*UP), stroke_color=WHITE)
        t2arcpoint1 = r_C2_O + 0.5*l2*UP
        t2arcpoint2 = r_C2_O + 0.5*l2*er2
        t2arc = ArcBetweenPoints(t2arcpoint2, t2arcpoint1)
        t2text = TexMobject(r'\theta_2')
        t2text.next_to(t2arc, UP)
        t2text.scale(1.25)
#        Line(np.array(r_C1_O), np.array(r_m1_O), stroke_color=GREY_1)
        self.play(Write(t2refline),
                  run_time=0.5)
        self.play(Write(t2arc),
                  run_time=0.5)
        self.play(Write(t2text))
        
        # x
        xrefline = DashedLine(np.array(6*LEFT+4*DOWN), np.array(6*LEFT+4*UP), stroke_color=WHITE)
        xcenterline = DashedLine(np.array(1*DOWN), np.array(5*DOWN), stroke_color=WHITE)
        xarrow = Arrow(np.array(6.23*LEFT+3.5*DOWN), np.array(0.26*RIGHT+3.5*DOWN))
        xtext = TexMobject(r'x')
        xtext.next_to(xarrow, UP)
        xtext.scale(1.25)
        self.play(Write(xrefline), Write(xcenterline), DrawBorderThenFill(xarrow),
                  run_time=0.5)
        self.play(Write(xtext),
                  run_time=0.5)
        self.wait(0.5)
        
        # Gravity
        g1arrow = Arrow(np.array(r_m1_O), np.array(r_m1_O+2*DOWN), stroke_color=NEON_GREEN)
        g1text = TexMobject(r'g')
        g1text.next_to(g1arrow, LEFT)
        g1text.scale(1.5)
        g2arrow = Arrow(np.array(r_m2_O), np.array(r_m2_O+2*DOWN), stroke_color=NEON_GREEN)
        g2text = TexMobject(r'g')
        g2text.next_to(g2arrow, RIGHT)
        g2text.scale(1.5)
        self.play(Write(g1arrow), Write(g1text),
                  Write(g2arrow), Write(g2text))
        self.wait(1)
        
        
        # Applied force F
        Farrow = Arrow(np.array(2.26*RIGHT+1*DOWN), np.array(5*RIGHT+1*DOWN), stroke_color=NEON_GREEN)
        Ftext = TexMobject(r'F')
        Ftext.next_to(Farrow, UP)
        Ftext.scale(1.5)
        self.play(Write(Farrow), Write(Ftext))
        self.wait(0.5)
        utext=TexMobject('u')
        utext.move_to(Ftext)
        utext.scale(1.5)
        self.play(Transform(Ftext, utext))
        self.wait(2)
        
        
        
        # Fade things to prepare to move stuff
        self.play(FadeOut(xrefline), FadeOut(xtext), FadeOut(xcenterline), FadeOut(xarrow),
                  FadeOut(t1refline), FadeOut(t1arc), FadeOut(t1text), 
                  FadeOut(t2refline), FadeOut(t2arc), FadeOut(t2text), 
                  FadeOut(Farrow), FadeOut(Ftext), 
                  FadeOut(g1arrow), FadeOut(g1text),
                  FadeOut(g2arrow), FadeOut(g2text),
                  FadeOut(l1text_1), FadeOut(l2text_1), 
                  FadeOut(m1text_1), FadeOut(m2text_1),
                  FadeOut(Mtext))

        
        
        # Demonstrate positional movement       
        self.play(cart.shift, 1*RIGHT,
                  MaintainPositionRelativeTo(p1, cart),
                  MaintainPositionRelativeTo(p2, cart))
        self.play(cart.shift, 2*LEFT,
                  MaintainPositionRelativeTo(p1, cart),
                  MaintainPositionRelativeTo(p2, cart))
        self.play(cart.shift, 1*RIGHT,
                  MaintainPositionRelativeTo(p1, cart),
                  MaintainPositionRelativeTo(p2, cart))
        
        self.wait(1)        
        self.play(Rotate(p1, t1/2, OUT, about_point=rod1.points[0]),
                  Rotate(p2, t1/2, OUT, about_point=rod2.points[0]))
        self.play(Rotate(p1, -t1/2, OUT, about_point=rod1.points[0]),
                  Rotate(p2, -t1/2, OUT, about_point=rod2.points[0]))
        
        self.wait(0.5)
        self.play(Rotate(p1, t1, OUT, about_point=rod1.points[0]),
                  Rotate(p2, t1, OUT, about_point=rod2.points[0]))
        self.wait(0.5)
        self.play(Rotate(p1, -2*t1, OUT, about_point=rod1.points[0]),
                  Rotate(p2, -2*t1, OUT, about_point=rod2.points[0]))
        self.wait(0.5)
        self.play(Rotate(p1, t1, OUT, about_point=rod1.points[0]),
                  Rotate(p2, t1, OUT, about_point=rod2.points[0]))
        self.play(FadeIn(l1text_1), FadeIn(l2text_1),
                  FadeIn(m1text_1), FadeIn(m2text_1))
        
        m1text_2 = TexMobject(r'm_1')
        m1text_2.move_to(m1text_1)
        m2text_2 = TexMobject(r'm_2')
        m2text_2.move_to(m2text_1)
        
        l1text_2 = TexMobject(r'l_1')
        l1text_2.move_to(l1text_1)
        l2text_2 = TexMobject(r'l_2')
        l2text_2.move_to(l2text_1)
        
        self.wait(0.5)
        self.play(Transform(m1text_1, m1text_2),
                  Transform(m2text_1, m2text_2),
                  Transform(l1text_1, l1text_2),
                  Transform(l2text_1, l2text_2))
        self.wait(0.5)
        self.play(m1.scale, 2,
                  m1text_1.shift, 0.2*UP)
        self.wait(0.5)
        self.play(rod2.stretch_about_point, 0.75, [1,0,0], rod2.points[0],
                  MaintainPositionRelativeTo(l2text_1, rod2),
                  m2.shift,-0.25*l2*er2,
                  MaintainPositionRelativeTo(m2text_1,m2))
        self.wait(0.5)
        self.play(FadeIn(xrefline), FadeIn(xtext), FadeIn(xcenterline), FadeIn(xarrow),
                  FadeIn(t1refline), FadeIn(t1arc), FadeIn(t1text), 
                  FadeIn(t2refline), FadeIn(t2arc), FadeIn(t2text), 
                  FadeIn(Farrow), FadeIn(Ftext),
                  FadeIn(Mtext))
        
        self.wait(3)
        
    def initialConditionPic(self):
        m1radius = 0.25
        m2radius = 0.25
        
        # Lengths, should also add angle
        l1 = 3.75
        l2 = 3.75
        
        # Angles (rad)
        t1 = PI/7
        t2 = PI/7
        
        # Unit vectors
        er1 = -np.sin(t1)*RIGHT+np.cos(t1)*UP
        et1 = -np.cos(t1)*RIGHT-np.sin(t1)*UP
        
        er2 = np.sin(t2)*RIGHT+np.cos(t2)*UP
        et2 = np.cos(t2)*RIGHT-np.sin(t2)*UP
        
        # Cart properties
        r_M_O = [0, -1] # Position of cart wrt the origin
        d_M = [0.75, 5] # Cart height, width
        
        # Wheel properties
        wheelradius = 0.4
        r_w1_O = [-3*wheelradius+r_M_O[0], -0.25*d_M[0]-wheelradius+r_M_O[1]]
        r_w2_O = [3*wheelradius+r_M_O[0], -0.25*d_M[0]-wheelradius+r_M_O[1]]
        
        # Ground properties
        r_G = [0,-0.5*d_M[0]-1.64*wheelradius+r_M_O[1],0]
        
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
        rod1 = Line(np.array(r_C1_O), np.array(r_m1_O), stroke_color=GREY_1)
        rod2 = Line(np.array(r_C2_O), np.array(r_m2_O), stroke_color=GREY_1)
        
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
        
        # text 
        Mtext = TexMobject(r'M')
        Mtext.move_to(M)
        
        m1text_1 = TexMobject(r'm')
        m1text_1.next_to(m1, UP)
        l1text_1 = TexMobject(r'l')
        l1text_1.next_to(rod1, LEFT)
        l1text_1.shift(0.75*RIGHT)
        
        m2text_1 = TexMobject(r'm')
        m2text_1.next_to(m2, UP)
        l2text_1 = TexMobject(r'l')
        l2text_1.next_to(rod2, RIGHT)
        l2text_1.shift(0.75*LEFT)
        
        m1text_2 = TexMobject(r'm_1')
        m1text_2.move_to(m1text_1)
        m2text_2 = TexMobject(r'm_2')
        m2text_2.move_to(m2text_1)
        
        l1text_2 = TexMobject(r'l_1')
        l1text_2.move_to(l1text_1)
        l2text_2 = TexMobject(r'l_2')
        l2text_2.move_to(l2text_1)
        
        vargroup = VGroup(m1text_2, m2text_2, l1text_2, l2text_2, Mtext)
        
        ground = Line(np.array(r_G+6*LEFT), np.array(r_G+6*RIGHT))
        
        # Create groups to move together
        cart = VGroup(M, w1, w2)
        
        p1 = VGroup(m1, rod1)
        p2 = VGroup(m2, rod2)
        sys = VGroup(cart, p1, p2)
        
        
        ##
        # Angles and position x
        ##
        # theta1
        t1refline = DashedLine(np.array(r_C1_O), np.array(r_C1_O+7*UP), stroke_color=WHITE)
        t1arcpoint1 = r_C1_O + 0.5*l1*UP
        t1arcpoint2 = r_C1_O + 0.5*l1*er1
        t1arc = ArcBetweenPoints(t1arcpoint1, t1arcpoint2)
        t1text = TexMobject(r'\theta_1')
        t1text.next_to(t1arc, UP)
        t1text.scale(1.25)
        t1group = VGroup(t1refline, t1arc,t1text)
        
        # theta2
        t2refline = DashedLine(np.array(r_C2_O), np.array(r_C2_O+7*UP), stroke_color=WHITE)
        t2arcpoint1 = r_C2_O + 0.5*l2*UP
        t2arcpoint2 = r_C2_O + 0.5*l2*er2
        t2arc = ArcBetweenPoints(t2arcpoint2, t2arcpoint1)
        t2text = TexMobject(r'\theta_2')
        t2text.next_to(t2arc, UP)
        t2text.scale(1.25)
        t2group = VGroup(t2refline, t2arc,t2text)
        
        # x
        xrefline = DashedLine(np.array(6*LEFT+4*DOWN), np.array(6*LEFT+4*UP), stroke_color=WHITE)
        xcenterline = DashedLine(np.array(1*DOWN), np.array(5*DOWN), stroke_color=WHITE)
        xarrow = Arrow(np.array(6.23*LEFT+3.5*DOWN), np.array(0.26*RIGHT+3.5*DOWN))
        xtext = TexMobject(r'x')
        xtext.next_to(xarrow, UP)
        xtext.scale(1.25)
        xgroup = VGroup(xrefline, xcenterline, xarrow, xtext)
        
        Farrow = Arrow(np.array(2.26*RIGHT+1*DOWN), np.array(5*RIGHT+1*DOWN), stroke_color=NEON_GREEN)
        Ftext = TexMobject(r'u')
        Ftext.next_to(Farrow, UP)
        Ftext.scale(1.5)
        Fgroup = VGroup(Farrow, Ftext)
        
        
        self.add(sys, ground, t1group, t2group, xgroup, Fgroup, vargroup)
        self.wait(0.5)
        
        l2factor = 0.9
        m2factor = 1
        m1density = 1/(PI*m2radius**2)
        r2 = np.sqrt(m2factor/(m1density*PI))
        m2scale = r2/m2radius
        t1f = 2.5 
        t2f = 2.5
        t1rot2ic = t1f*PI/180 - t1
        t2rot2ic = t2 - t2f*PI/180
        
        l2reltext = TexMobject(str(round(l2factor,1)), r'l_1')
        l2reltext.move_to(l2text_2)
        l2reltext.shift(0.3*RIGHT)
        if m2factor == 1:
            m2reltext = TexMobject(r'm_1')
        else:
            m2reltext = TexMobject(str(round(m2factor,1)), r'm_1')
        
        # Change the length
        self.play(ReplacementTransform(l2text_2, l2reltext))
        self.play(rod2.stretch_about_point, l2factor, [1,0,0], rod2.points[0],
                  t2arc.stretch_about_point, l2factor, [1,0,0], rod2.points[0],
                  t2text.stretch_about_point, l2factor, [1,0,0], rod2.points[0],
#                  l2reltext.stretch_about_point, l2factor, [1,0,0], rod2.points[0],
                  MaintainPositionRelativeTo(l2reltext, rod2),
                  m2.shift,-(1-l2factor)*l2*er2,
                  MaintainPositionRelativeTo(m2text_2,m2))
        self.wait(0.5)
        
        # Change mass
        m2reltext.move_to(m2text_2)
        self.play(ReplacementTransform(m2text_2, m2reltext))
#        self.play(m2.scale, 1.5)
        
        # Change angles
        t1ic = TexMobject(r'\theta_1(0)', r'=', str(t1f), r'^{\circ}')
        t2ic = TexMobject(r'\theta_2(0)', r'=', str(t2f), r'^{\circ}')
        t1ic.next_to(xrefline, RIGHT)
        t1ic.set_y(1)
        t2ic.next_to(t1ic, DOWN)
        t2ic.align_to(t1ic, LEFT)
        
        self.play(ReplacementTransform(t1text, t1ic), ReplacementTransform(t2text, t2ic),
                  FadeOut(t1arc), FadeOut(t2arc))
        self.wait(0.5)
        self.play(Rotate(p1, t1rot2ic, OUT, about_point=rod1.points[0]),
                  Rotate(p2, t2rot2ic, OUT, about_point=rod2.points[0]),
                  MaintainPositionRelativeTo(m2reltext, m2),
                  MaintainPositionRelativeTo(m1text_2, m1),
                  MaintainPositionRelativeTo(l1text_2, rod1),
                  MaintainPositionRelativeTo(l2reltext, rod2),)
        self.wait(1)
        
        l1val = TexMobject(r'l_1 = 2 m')
        m1val = TexMobject(r'm_1 = 0.1 kg')
        Mval = TexMobject(r'M = 1 kg')
        Mval.next_to(t1ic, UP)
        l1val.next_to(Mval, UP)
        m1val.next_to(l1val, UP)
        Mval.align_to(t1ic, LEFT)
        m1val.align_to(t1ic, LEFT)
        l1val.align_to(t1ic, LEFT)
        self.play(TransformFromCopy(Mtext, Mval),
                  TransformFromCopy(m1text_2, m1val),
                  TransformFromCopy(l1text_2, l1val))
        
        self.wait(2)
        
        # Next condition
        
        
        
    def getPos(self, obj):
        objpos = [obj.get_x(), obj.get_y(), obj.get_z()]
        return objpos
            
    def relVec(self, obj, pos):
        # generate the vector to move an object to a new position
        return pos+self.getPos(obj)
    
    
#    class Pendulum(Circle, Line, radius, length, mass, angle, position):
#        CONFIG = {
#                "radius": radius,
#                "fill_color": GREY_1,
#                "stroke_color": NEON_GREEN,
#                "fill_opacity": 0.5,
#                }
#        def __init__(self, **kwargs):
#            Circle.__init__(**kwargs)