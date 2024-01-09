from sympy import symbols, Matrix, sqrt, sin, cos, log, tan, exp, simplify, collect, conjugate, expand
from sympy.abc import f, F
import numpy as np

# symbols
etaM, etam, fiM, fim, eta0, fi0 = symbols('etaM etam fiM fim eta0 fi0')
y1, y2, y3, alfa, beta, sig = symbols('y1 y2 y3 alfa beta sig')
h0, H0, A0, G0, HM, Hm, GM, Gm = symbols('h0 H0 A0 G0 HM Hm GM Gm')
Tb53 = symbols(' Tb53')
T6, T5, t, b, T, Tb23, Bp = symbols('T6 T5 t b T Tb23 Bp')
tbar, T6bar, T5bar, bbar, Tbar, Tb23bar, Bbar = symbols('tbar T6bar T5bar bbar Tbar Tb23bar Bbar')
Tb53bar, u, d = symbols('Tb53bar u d')
ubar, dbar, c, s, cbar, sbar = symbols('ubar dbar c s cbar sbar')
pr, yb, yB, yd, ys, yu, yc = symbols('pr yb yB yd ys yu yc')


def scalar_lagrangian():
    '''
    Method to calculate the scalar lagrangian for the Bestets Little Higgs Model
    '''
    # Constants
    a = 1j/sqrt(2)
    co = 1j * (1/2)
    r = 1/sqrt(2)
    v = 246
    x = v/f
    y = v**2 / f**2
    y12 = 1/sqrt(y1**2 + y2**2)
    y13 =1/sqrt(y1**2 + y3**2)
    y1m2 = 1/(y1**2 + y2**2)
    y1m3 = 1/(y1**2 + y3**2)
    y2My3 = 1/(y2**2 - y3**2)
    y1d3 = 2*y1**2 - y3**2
    A = 1/sqrt(y1**2 + y2**2)
    B = 1/sqrt(y1**2 + y3**2)
    cb = cos(beta)
    sb = sin(beta)
    ta = tan(alfa)
    ca = cos(alfa)
    sa = sin(alfa)
    f2 = 1/sqrt(f**2 + F**2)

    fi1 = f2*r*F*(fiM + fim)
    fi2 = -1j*r*f2*F*(fiM - fim)
    fi3 = fi0*f2*F
    fi2p = 1j*r*f2*F*(fiM - fim)
    eta1 = r*(etaM + etam)
    eta2 = 1j*r*(etaM - etam)
    eta2p = -1j*r*(etaM + etam)
    eta3 = eta0

    h11 = ca*h0-sa*H0+v*sb
    h12 = -cb*A0+sb*G0
    h13 = -r*(cb*(HM + Hm)-sb*(GM + Gm))
    h14 = 1j*r*(cb*(HM - Hm) - sb*(GM - Gm))
    h21 = sa*h0 + ca*H0 + v*cb
    h22 = sb*A0 + cb*G0
    h23 = r*(sb*(HM +Hm) + cb*(GM + Gm))
    h24 = -1j*r*(sb*(HM - Hm) + cb*(GM - Gm))

    hT11 = ca*h0 -sa*H0 +v*sb
    hT12 = -cb*A0 + sb*G0
    hT13 = -r*(cb*(HM + Hm)-sb*(GM + Gm))
    hT14 = 1j*r*(cb*(HM - Hm) - sb*(GM - Gm))
    hT21 = sa*h0 + ca*H0 + v*cb
    hT22 = sb*A0 + cb*G0
    hT23 = r*(sb*(HM + Hm) + cb*(GM + Gm))
    hT24 = -1j*r*(sb*(HM - Hm) + cb*(GM - Gm))

    TR2 = np.array([[0,0,-co,0,0,0],[0,0,0,-co,0,0],[co,0,0,0,0,0],[0,co,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])

    PI1 = np.array([[0,0,0,co*(-eta1+fi1),0,0],[0,0,co*(eta1+fi1),0,0,0],
                    [0,-co *(eta1+fi1),0,0,0,0],[co*(eta1-fi1),0,0,0,0,0],
                    [0,0,0,0,0,0],[0,0,0,0,0,0]])
    PI2 = np.array([[0,0,co*(fi2-eta2),0,0,0],[0,0,0,-co*(eta2+fi2),0,0],
                [co*(eta2-fi2),0,0,0,0,0],[0,co*(eta2+fi2),0,0,0,0],
                [0,0,0,0,0,0],[0,0,0,0,0,0]])
    PI3 = np.array([[0,co*(fi3-eta3),0,0,0,0],[co*(eta3-fi3),0,0,0,0,0],
                    [0,0,0,co*(eta3+fi3),0,0],[0,0,-co*(eta3+fi3),0,0,0],
                    [0,0,0,0,0,a*sig],[0,0,0,0,-a*sig,0]])

    PI = PI1 + PI2 + PI3
    vh1 = np.array([0,0,0,0,a*h11,a*h21])
    vh2 = np.array([0,0,0,0,a*h12,a*h22])
    vh3 = np.array([0,0,0,0,a*h13,a*h23])
    vh4 = np.array([0,0,0,0,a*h14,a*h24])
    vh5 = np.array([-a*hT11,-a*hT12,-a*hT13,-a*hT14,0,0])
    vh6 = np.array([-a*hT21,-a*hT22,-a*hT23,-a*hT24,0,0])
    PIh = np.array([vh1,vh2,vh3,vh4,vh5,vh6])
    SumaP = PI + PIh

    Uno6 = np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
    P5 = np.array([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,1,0],[0,0,0,0,0,0]])
    P6 = np.array([[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,1]])
    S = np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,-1,0],[0,0,0,0,0,-1]])

    aux1 = PI + PIh
    aux2 = Uno6+(2j/f)*aux1
    aux3 = aux1 @ aux1
    SigmaP = aux2 - (2/f**2)*aux3

    PiT = np.transpose(-PI)

    ht11 = ca*h0 - sa*H0 + v*sb
    ht12 = -cb*A0 + sb*G0
    ht13 = -r*(cb*(HM + Hm) - sb*(GM + Gm))
    ht14 = -1j*r*(cb*(Hm - HM) - sb*(Gm - GM))
    ht21 = sa*h0 + ca*H0 + v*cb
    ht22 = sb*A0 + cb*G0
    ht23 = r*(sb*(HM + Hm) + cb*(GM + Gm))
    ht24 = 1j*r*(sb*(Hm - HM)+cb*(Gm - GM))

    htT11 = ca*h0 - sa*H0 + v*sb
    htT12 = -cb*A0 + sb*G0
    htT13 = -r*(cb*(HM + Hm) -sb*(GM + Gm))
    htT14 = -1j*r*(cb*(Hm - HM) -sb*(Gm - GM))
    htT21 = sa*h0 + ca*H0 + v*cb
    htT22 = sb*A0 + cb*G0
    htT23 = r*(sb*(HM + Hm) + cb*(GM + Gm))
    htT24 = 1j*r*(sb*(Hm -HM) + cb*(Gm - GM))

    vht1 = np.array([0,0,0,0,a*htT11,a*htT21])
    vht2 = np.array([0,0,0,0,a*htT12,a*htT22])
    vht3 = np.array([0,0,0,0,a*htT13,a*htT23])
    vht4 = np.array([0,0,0,0,a*htT14,a*htT24])
    vht5 = np.array([-a*ht11,-a*ht12,-a*ht13,-a*ht14,0,0])
    vht6 = np.array([-a*ht21,-a*ht22,-a*ht23,-a*ht24,0,0])
    PihT = np.array([vht1,vht2,vht3,vht4,vht5,vht6])

    aux11 = PihT + PiT
    aux22 = Uno6+(2j/f)*aux11
    aux33 = aux11 @ aux11
    SigmaPT = aux22 - (2/f**2)*aux33

    # Parametrization of the fundamental fields Q withind 6D
    q31 = t-2*x*cb*y2*y12*T6 - x*sb*y2*(2*y1**2 - y3**2)*y12*y1m3*T5
    Qu = T-2*x*cb*y1*y12*T6-x*sb*y1*(2*y2**2+y3**2)*y12*y2My3*T5
    q32 = b
    Qb23 = Tb23+x*sb*T5
    Q5 = T5-x*sb*Tb23+x*sb*y2*(2*y1**2-y3**2)*y12*y1m3*t+x*sb*y1*(2*y2**2+y3**2)*y12*y2My3*T
    Q6 = T6+2*x*cb*y1*y12*T+2*x*cb*y2*y12*t
    Qd = Bp
    Qb53 = Tb53
    Qa1 = y12*(y1*Qu+y2*q31)
    Qpa1 = y12*(y2*Qu-y1*q31)
    Qa2 = y12*(y1*Qd+y2*q32)
    Qpa2 = y12*(y2*Qd-y1*q32)
    Qb1 = Qb53
    Qb2 = Qb23

    # transposes
    q31bar = tbar-2*x*cb*y2*y12*T6bar-x*sb*y2*(2*y1**2-y3**2)*y12*y1m3*T5bar
    q32bar = bbar
    Qubar = Tbar-2*x*cb*y1*y12*T6bar-x*sb*y1*(2*y2**2+y3**2)*y12*y2My3*T5bar
    Qb23bar = Tb23bar+x*sb*T5bar
    Q5bar = T5bar-x*sb*Tb23bar+x*sb*y2*(2*y1**2-y3**2)*y12*y1m3*tbar+x*sb*y1*(2*y2**2+y3**2)*y12*y2My3*Tbar
    Q6bar = T6bar+2*x*cb*y1*y12*Tbar+2*x*cb*y2*y12*tbar
    Qdbar = Bbar
    Qb53bar = Tb53bar
    Qa1bar = y12*(y1*Qubar+y2*q31bar)
    Qpa1bar = y12*(y2*Qubar-y1*q31bar)
    Qa2bar = y12*(y1*Qdbar+y2*q32bar)
    Qpa2bar = y12*(y2*Qdbar-y1*q32bar)
    Qb1bar = Qb53bar
    Qb2bar = Qb23bar

    # doublets (leptones absents)
    q1 = np.transpose([-r*u,1j*r*u,r*d,1j*r*d,0,0])
    q1bar = [-r*ubar,1j*r*ubar,r*dbar,1j*r*dbar,0,0]
    q2 = np.transpose([-r*c,1j*r*c,r*s,1j*r*s,0,0])
    q2bar = [-r*cbar,1j*r*cbar,r*sbar,1j*r*sbar,0,0]
    q3 = np.transpose([-r*t,1j*r*t,r*b,1j*r*b,0,0])
    q3bar = [-r*tbar,1j*r*tbar,r*bbar,1j*r*bbar,0,0]

    # singlets (leptons absents)
    u1 = np.transpose([0,0,0,0,u,0])
    u1bar = [0,0,0,0,ubar,0]
    u2 = np.transpose([0,0,0,0,c,0])
    u2bar = [0,0,0,0,cbar,0]
    d1 = np.transpose([0,0,0,0,d,0])
    d1bar = [0,0,0,0,dbar,0]
    d2 = np.transpose([0,0,0,0,s,0])
    d2bar = [0,0,0,0,sbar,0]
    d3 = np.transpose([0,0,0,0,b,0])
    d3bar = [0,0,0,0,bbar,0]
    d4 = np.transpose([0,0,0,0,Bp,0])
    d4bar = [0,0,0,0,Bbar,0]

    # Fundamental Fields Q within 6D
    Q = np.transpose([-r*(Qa1-Qb2),1j*r*(Qa1-Qb2),r*(Qa2-Qb1),1j*r*(Qa2+Qb1),r*sqrt(2)*Q5,r*sqrt(2)*Q6])
    Qbar = [r*(-Qa1bar-Qb2bar),1j*r*(Qa1bar-Qb2bar),r*(Qa2bar-Qb1bar),1j*r*(Qa2bar+Qb1bar),r*sqrt(2)*Q5bar,r*sqrt(2)*Q6bar]
    Qp = np.transpose([-r*Qpa1,1j*r*Qpa1,r*Qpa2,1j*r*Qpa2,0,0])
    Qpbar = [-r*Qpa1bar,1j*r*Qpa1bar,r*Qpa2bar,1j*r*Qpa2bar,0,0]

    # parametrizaciones para los campos fundamentales U
    Ua = T-x*cb*T6+x*(sb*y3*(2*y1**2-y2**2))*(1/(2*y1m2**2*y1m3))*t+x*(sb*y1*(y2**2+2*y3**2)*y13)*(1/(y3**2-y2**2))*T5
    u3 = t-x*2*sb*y3*y13*Tb23-x*sb*y3*(2*y1**2-y2**2)*y13*y1m2*T
    Uu = T-x*cb*T6+x*sb*y3*(2*y1**2-y2**2)*2*y1m2**2*(1/(-y2**2+y3**2))*t+x*sb*y1*(y2**2+2*y3**2)*y13*(1/(-y2**2+y3**2))*T5
    Ub23 = Tb23+2*x*sb*y1*y13*T5+2*x*sb*y3*y13*t
    Uq5 = T5+2*x*sb*y1*y13*Tb23+2*x*sb*y3*y13*t
    U6 = T6+x*cb*T
    Ud = Bp
    Ub53 = Tb53
    Ua1 = -Ud
    Ua2 = Ua
    Ub1 = Ub23
    Ub2 = -Ub53
    U5 = y13*(y1*Uq5+y3*u3)
    Up5 = y13*(y3*Uq5-y1*u3)

    # Transposes
    Uabar = Tbar-x*cb*T6bar+x*(sb*y3*(2*y1**2-y2**2))*(1/(2*y1m2**2*y1m3))*tbar+x*(sb*y1*(y2**2+2*y3**2)*y13)*(1/(y3**2-y2**2))*T5bar
    u3bar = tbar-x*2*sb*y3*y13*Tb23bar-x*sb*y3*(2*y1**2-y2**2)*y13*y1m2*Tbar
    Uubar = Tbar-x*cb*T6bar+x*sb*y3*(2*y1**2-y2**2)*2*y1m2**2*(1/(-y2**2+y3**2))*tbar+x*sb*y1*(y2**2+2*y3**2)*y13*(1/(-y2**2+y3**2))*T5bar
    Ub23bar = Tb23bar+2*x*sb*y1*y13*T5bar+2*x*sb*y3*y13*tbar
    Uq5bar = T5bar+2*x*sb*y1*y13*Tb23bar+2*x*sb*y3*y13*tbar
    U6bar = T6bar+x*cb*Tbar
    Udbar = Bbar
    Ub53bar = Tb53bar
    Ua1bar = -Udbar
    Ua2bar = Uabar
    Ub1bar = Ub23bar
    Ub2bar = -Ub53bar
    U5bar = y13*(y1*Uq5bar+y3*u3bar)
    Up5bar = y13*(y3*Uq5bar-y1*u3bar)

    # Fundamental Fields U within 6D
    Uc = np.transpose([-r*Ub1-r*Ua2,1j*r*(Ub1-Ua2),r*(Ub2-Ua1),1j*r*(Ub2+Ua1),r*sqrt(2)*U5,r*sqrt(2)*U6])
    Ucbar = [-r*Ub1bar-r*Ua2bar,1j*r*(Ub1bar-Ua2bar),r*Ub2bar-r*Ua1bar,1j*r*(Ub2bar+Ua1bar),r*sqrt(2)*U5bar,r*sqrt(2)*U6bar]
    Upc5 = np.transpose([0,0,0,0,Up5,0])
    Upc5bar = [0,0,0,0,Up5bar,0]

    # Right Lagrangians (partials)
    L1r = pr*y1*f*(((Qbar @ S) @ (SigmaP @ S)) @ Uc)
    L2r = pr*y2*f*((Qpbar @ SigmaP) @ Uc)
    L3r = pr*y3*f*((Qbar @ SigmaP) @ Upc5)
    L4r = -2j*pr*yb*f*((q3bar @ TR2) @ (SigmaP @ d3))
    L12r = -2j*pr*yB*f*((q1bar @ TR2) @ (SigmaP @ d4))
    L13r = -2j*pr*yB*f*((q2bar @ TR2) @ (SigmaP @ d4))
    L5r = -2j*pr*yd*f*((q1bar @ TR2) @ (SigmaP @ d1))
    L6r = -2j*pr*ys*f*((q2bar @ TR2) @ (SigmaP @ d2))
    L10r = pr*yu*f*(q1bar @ (SigmaP @ u1))
    L11r = pr*yc*f*(q2bar @ (SigmaP @ u2))

    # Lagrangiano derecho
    Lr1 = L1r+L2r+L3r+L4r+L5r+L6r+L10r+L11r+L12r+L13r

    pl = symbols('pl')

    # Left Lagrangians (partials)
    L1d = pl*y1*f*((Ucbar @ S) @ (SigmaPT @ S) @ Q)
    L2d = pl*y2*f*((Ucbar @ SigmaPT) @ Qp)
    L3d = pl*y3*f*((Upc5bar @ SigmaPT) @ Q)
    L4d = 2j*pl*yb*f*((d3bar @ SigmaPT) @ (TR2 @ q3))
    L5d = 2j*pl*yd*f*((d1bar @ SigmaPT) @ (TR2 @ q1))
    L6d = 2j*pl*ys*f*((d2bar @ SigmaPT) @ (TR2 @ q2))
    L10d = pl*yu*f*((u1bar @ SigmaPT) @ q1)
    L11d = pl*yc*f*((u2bar @ SigmaPT) @ q2)
    L12d = 2j*pl*yB*f*((d4bar @ SigmaPT) @ (TR2 @ q1))
    L13d = 2j*pl*yB*f*((d4bar @ SigmaPT) @ (TR2 @ q2))    
    Lt1 = L1d + L2d +L3d + L4d + L5d + L6d + L10d + L11d + L12d + L13d
    
    # Total scalar lagrangian
    l_scalar = Lr1 + Lt1
    return l_scalar

def scalar_vertex(lag, qbar, field, q): #, beta_value=0, alfa_value=np.pi/2,
                 #f_value=1000, y1_value=0.9, y2_value=0.7, y3_value=0.33):
    '''
    Method to obtain the vertex
    @Parameters: vertex format = qbar, field, q
    '''
    q1 = str(q)
    field1 = str(field)
    qbar1 = str(qbar)

    # aqui ya queda el lagrangiano sin qbar, field y q
    quark = ["t","T","T5","T6","Tb23","Tb53","Bp"]
    quarkbar = ["tbar","Tbar","T5bar","T6bar","Tb23bar","Tb53bar","Bbar"]
    fields = ["A0","H0","h0","sig","HM","Hm","fi0","fiM","fim","eta0","etaM","etam"]
    value_q = [0, 0, 0, 0, 0, 0, 0]
    value_qbar = [0, 0, 0, 0, 0, 0, 0]
    value_field = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
    i = 0
    for fila in quark:
        i = i + 1
        if q1 == quark[i - 1]:
            indice_q = i - 1
            value_q[indice_q] = q
            # print("\n indice: ", i-1, value_q)
        else:
            continue

    j = 0
    for fila in quarkbar:
        j = j + 1
        if qbar1 == quarkbar[j - 1]:
            indice_qbar = j - 1
            value_qbar[indice_qbar] = qbar
            # print("\n indice: ", j - 1, value_qbar)
        else:
            continue
    k = 0
    for fila in fields:
        k = k + 1
        if field1 == fields[k - 1]:
            indice_fields = k - 1
            value_field[indice_fields] = field
            # print("\n indice: ", k - 1, value_field)
        else:
            continue
    # Evaluates the lagrangian
    Lagrangiano_evaluado = lag.subs({d: 0, dbar: 0, u: 0, ubar: 0, s:0, sbar: 0, c: 0, cbar: 0, b: 0,
                                     bbar: 0, t: value_q[0], T: value_q[1], T5: value_q[2], T6: value_q[3],
                                     Tb23: value_q[4], Tb53: value_q[5], Bp: value_q[6], tbar: value_qbar[0],
                                     Tbar: value_qbar[1], T5bar: value_qbar[2], T6bar: value_qbar[3],
                                     Tb23bar: value_qbar[4], Tb53bar: value_qbar[5], Bbar: value_qbar[6],
                                     A0: value_field[0], H0: value_field[1], h0: value_field[2],
                                     sig: value_field[3], HM: value_field[4], Hm: value_field[5],
                                     fi0: value_field[6], fiM: value_field[7], fim: value_field[8],
                                     eta0: value_field[9], etaM: value_field[10], etam: value_field[11],
                                     GM: 0, Gm: 0, G0: 0})
    # print(Lagrangiano_evaluado) # shows the lagrangian only with the contributions: qbar,field,q
    Lag_expan = expand(Lagrangiano_evaluado) # we need to expand at every step to avoid errors
    expr_q = collect(Lag_expan, q).coeff(q, 1)
    expr_q_expan = expand(expr_q)
    expr_qbar = collect(expr_q_expan, qbar).coeff(qbar, 1)
    expr_qbar_expan = expand(expr_qbar)
    expr_field = collect(expr_qbar_expan, field).coeff(field, 1)
    # print(expr_field) # prints the vertex
    return expr_field


def scalar_vertex_value(expr_field, beta_value=0, alfa_value=np.pi/2,
                        f_value=1000, y1_value=0.9, y2_value=0.7, y3_value=0.33):
    '''
    Method to obtain the numerical value of the vertex
    '''
    vertex_value = expr_field.subs({beta: beta_value, alfa: alfa_value, f: f_value, y1: y1_value, y2: y2_value, y3: y3_value, sqrt(2): 1.41421})
    return vertex_value # return type=sympycore object
