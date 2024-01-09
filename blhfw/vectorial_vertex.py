import numpy as np
from sympy import symbols, cos, sin, sqrt, collect, expand

beta, alfa, f, y1, y2, y3, gB = symbols('beta alfa f y1 y2 y3 gB')
t, T, T5, T6, Tb23, Tb53, Bp, b = symbols('t T T5 T6 Tb23 Tb53 Bp b')
tbar, Tbar, T5bar, T6bar, Tb23bar, Tb53bar, Bbar, bbar = symbols('tbar Tbar T5bar T6bar Tb23bar Tb53bar Bbar bbar')
Z, Zp, WM, Wm, WpM, Wpm, gamma = symbols('Z Zp WM Wm WpM Wpm gamma')
u, c, s, d = symbols('u c s d')
ubar, cbar, sbar, dbar = symbols('ubar cbar sbar dbar')
cp = 1
gp = 0.3528
g = 0.6528
gA = gB


def vector_lagrangian():
    '''
    Method to calculate the scalar lagrangian for the Bestets Little Higgs Model
    '''
    # CONSTANTS
    #beta = 1.35
    #alfa = 0.15
    cp = 1
    r = 1/np.sqrt(2)
    v = 246
    gp = 0.3528
    g = 0.6528
    gA = gB
    tw = np.arctan(gp/g)
    sw = np.sin(tw)
    cw = np.cos(tw)
    cb = cos(beta)
    sb = sin(beta)
    ca = cos(alfa)
    sa = sin(alfa)
    cg = g/gA
    sg = g/gB
    co = 1j/2.0

    # TERMINOS REPETIDOS
    x = v/f
    y = v**2/f**2
    y12 = 1/sqrt(y1**2 + y2**2)
    y13 = 1/sqrt(y1**2 + y3**2)
    y1m2 = 1/(y1**2 + y2**2)
    y1m3 = 1/(y1**2 + y3**2)
    y2My3 = 1/(y2**2-y3**2)
    y1d3 = 2*y1**2 - y3**2
    A = 1/sqrt(y1**2 + y2**2)
    B = 1/sqrt(y1**2 + y3**2) # es el primer B

    # PARAMETRIZACION DE LOS Qs
    q31 = t - 2*x*cb*y2*y12*T6 - x*sb*y2*(2*y1**2 - y3**2)*y12*y1m3*T5
    q32 = b
    Qu = T - 2*x*cb*y1*y12*T6 - x*sb*y1*(2*y2**2 + y3**2)*y12*y2My3*T5
    Qb23 = Tb23 + x*sb*T5
    Q5 = T5 - x*sb*Tb23 + x*sb*y2*(2*y1**2 - y3**2)*y12*y1m3*t + x*sb*y1*(2*y2**2 + y3**2)*y12*y2My3*T
    Q6 = T6 + 2*x*cb*y1*y12*T + 2*x*cb*y2*y12*t
    Qd = Bp
    Qb53 = Tb53

    q31bar = tbar - 2*x*cb*y2*y12*T6bar - x*sb*y2*(2*y1**2 - y3**2)*y12*y1m3*T5bar
    q32bar = bbar
    Qubar = Tbar - 2*x*cb*y1*y12*T6bar - x*sb*y1*(2*y2**2 + y3**2)*y12*y2My3*T5bar
    Qb23bar = Tb23bar + x*sb*T5bar
    Q5bar = T5bar - x*sb*Tb23bar + x*sb*y2*(2*y1**2 - y3**2)*y12*y1m3*tbar + x*sb*y1*(2*y2**2 + y3**2)*y12*y2My3*Tbar
    Q6bar = T6bar + 2*x*cb*y1*y12*Tbar + 2*x*cb*y2*y12*tbar
    Qdbar = Bbar
    Qb53bar = Tb53bar

    Qa1 = y12*(y1*Qu + y2*q31)
    Qpa1 = y12*(y2*Qu - y1*q31)
    Qa2 = y12*(y1*Qd + y2*q32)
    Qpa2 = y12*(y2*Qd - y1*q32)
    Qb1 = Qb53
    Qb2 = Qb23

    Qa1bar = y12*(y1*Qubar + y2*q31bar)
    Qpa1bar = y12*(y2*Qubar - y1*q31bar)
    Qa2bar = y12*(y1*Qdbar + y2*q32bar)
    Qpa2bar = y12*(y2*Qdbar - y1*q32bar)
    Qb1bar = Qb53bar
    Qb2bar = Qb23bar

    a11 = r*(cg*(Wm + WM) + sg*(Wpm + WpM))
    a12 = 1j*r*(cg*(-Wm + WM) + sg*(-Wpm + WpM))
    a21 = r*(sg*(Wm + WM) - sg*(Wpm + WpM))
    a22 = 1j*r*(sg*(-Wm + WM) - cg*(-Wpm + WpM))
    a13 = cg*sw*gamma + cg*cw*cp*Z + sg*cp*Zp
    a23 = sg*sw*gamma + sg*cw*cp*Z - cg*cp*Zp
    b3 = cw*gamma -sw*cp*Z
    
    Uno6 = np.array([[1,0,0,0,0,0],[0,1,0,0,0,0],[0,0,1,0,0,0],[0,0,0,1,0,0],[0,0,0,0,1,0],[0,0,0,0,0,1]])
    TL3 = np.array([[0,co,0,0,0,0],[-co,0,0,0,0,0],[0,0,0,co,0,0],[0,0,-co,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
    TR3 = np.array([[0,-co,0,0,0,0],[co,0,0,0,0,0],[0,0,0,co,0,0],[0,0,-co,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
    TL1 = np.array([[0,0,0,co,0,0],[0,0,co,0,0,0],[0,-co,0,0,0,0],[-co,0,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])
    TL2 = np.array([[0,0,co,0,0,0],[0,0,0,-co,0,0],[-co,0,0,0,0,0],[0,co,0,0,0,0],[0,0,0,0,0,0],[0,0,0,0,0,0]])

    Q = np.transpose([-r*(Qa1+Qb2),1j*r*(Qa1-Qb2),r*(Qa2-Qb1),1j*r*(Qa2+Qb1),Q5,Q6])
    Qbar = np.array([-r*(Qa1bar+Qb2bar),-1j*r*(Qa1bar-Qb2bar),r*(Qa2bar-Qb1bar),-1j*r*(Qa2bar+Qb1bar),Q5bar,Q6bar])
    Qp = np.transpose([-r*Qpa1,1j*r*Qpa1,r*Qpa2,1j*r*Qpa2,0,0])
    Qpbar = np.array([-r*Qpa1bar,-1j*r*Qpa1bar,r*Qpa2bar,-1j*r*Qpa2bar,0,0])
    
    q1 = np.transpose([-r*u,1j*r*u,r*d,1j*r*d,0,0])
    q1bar = np.array([-r*ubar,-1j*r*ubar,r*dbar,-1j*r*dbar,0,0])
    q2 = np.transpose([-r*c,1j*r*c,r*s,1j*r*s,0,0])
    q2bar = np.array([-r*cbar,-1j*r*cbar,r*sbar,-1j*r*sbar,0,0])
    q3 = np.transpose([0,0,r*Bp,1j*r*Bp,0,0])
    q3bar = np.array([0,0,r*Bbar,-1j*r*Bbar,0,0])

    u1 = np.transpose([0,0,0,0,u,0])
    u1bar = np.array([0,0,0,0,ubar,0])
    u2 = np.transpose([0,0,0,0,c,0])
    u2bar = np.array([0,0,0,0,cbar,0])

    d1 = np.transpose([0,0,0,0,d,0])
    d1bar = np.array([0,0,0,0,dbar,0])
    d2 = np.transpose([0,0,0,0,s,0])
    d2bar = np.array([0,0,0,0,sbar,0])

    d3 = np.transpose([0,0,0,0,b,0])
    d3bar = np.array([0,0,0,0,bbar,0])

    Uno6p = (2/3)*Uno6
    DQ = 1j*gA*a11*(TL1 @ Q) + 1j*gA*a12*(TL2 @ Q) + 1j*gA*a13*(TL3 @ Q) + 1j*gp*b3*((TR3 @ Q) + (Uno6p @ Q))    
    DQp = 1j*gA*a11*(TL1 @ Qp) + 1j*gA*a12*(TL2 @ Qp) + 1j*gA*a13*(TL3 @ Qp) + 1j*gp*b3*((TR3 @ Qp) + (Uno6p @ Qp))
    Dq1 = 1j*gA*a11*(TL1 @ q1) + 1j*gA*a12*(TL2 @ q1) + 1j*gA*a13*(TL3 @ q1) + 1j*gp*b3*((TR3 @ q1) + (Uno6p @ q1))
    Dq2 = 1j*gA*a11*(TL1 @ q2) + 1j*gA*a12*(TL2 @ q2) + 1j*gA*a13*(TL3 @ q2) + 1j*gp*b3*((TR3 @ q2) + (Uno6p @ q2))    
    Du1 = 1j*gp*b3*(-2/3)*(Uno6 @ u1)
    Du2 = 1j*gp*b3*(-2/3)*(Uno6 @ u2)
    Dd1 = 1j*gp*b3*(1/3)*(Uno6 @ d1)
    Dd2 = 1j*gp*b3*(1/3)*(Uno6 @ d2)

    # LAGRANGIANOS DERECHOS
    L1d = 1j*(Qbar @ DQ)
    L2d = 1j*(Qpbar @ DQp)
    L9a = 1j*(q1bar @ Dq1)
    L10a = 1j*(q2bar @ Dq2)
    L101a = 1j*(q3bar @ Dq2)
    L91a = 1j*(q3bar @ Dq1)
    L11a = 1j*(u1bar @ Du1)
    L12a = 1j*(u2bar @ Du2)
    L13a = 1j*(d1bar @ Dd1)
    L14a = 1j*(d2bar @ Dd2)
    Ld = -L1d - L2d - L9a - L10a - L101a - L91a + L11a + L12a + L13a + L14a
    
    Ua = T - x*cb*T6 + x*(sb*y3*(2*y1**2 - y2**2))*(1/(y1m3*2*y1m2**2))*t + x*(sb*y1*(y2**2 + 2*y3**2)*y13)*(1/(y3**2 - y2**2))*T5
    u3 = t - x*2*sb*y3*y13*Tb23 - x*sb*y3*(2*y1**2 - y2**2)*y13*y1m2*T
    Uu = T - x*cb*T6 + x*sb*y3*(2*y1**2 - y2**2)*2*y1m2**2*(1/(-y2**2 + y3**2))*t + x*sb*y1*(y2**2 + 2*y3**2)*y13*(1/(-y2**2+y3**2))*T5
    Ub23 = Tb23 + 2*x*sb*y1*y13*T5 + 2*x*sb*y3*y13*t
    Uq5 = T5 + 2*x*sb*y1*y13*Tb23 + 2*x*sb*y3*y13*t
    U6 = T6 + x*cb*T
    Ud = Bp
    Ub53 = Tb53
    Uabar = Tbar - x*cb*T6bar + x*(sb*y3*(2*y1**2 - y2**2))*(1/(y1m3*2*y1m2**2))*tbar + x*(sb*y1*(y2**2 + 2*y3**2)*y13)*(1/(y3**2 - y2**2))*T5bar
    u3bar = tbar - x*2*sb*y3*y13*Tb23bar - x*sb*y3*(2*y1**2 - y2**2)*y13*y1m2*Tbar
    Uubar = Tbar - x*cb*T6bar + x*sb*y3*(2*y1**2 - y2**2)*2*y1m2**2*(1/(-y2**2 + y3**2))*tbar + x*sb*y1*(y2**2 + 2*y3**2)*y13*(1/(-y2**2+y3**2))*T5bar
    Ub23bar = Tb23bar + 2*x*sb*y1*y13*T5bar + 2*x*sb*y3*y13*tbar
    Uq5bar = T5bar + 2*x*sb*y1*y13*Tb23bar + 2*x*sb*y3*y13*tbar
    U6bar = T6bar + x*cb*Tbar
    Udbar = Bbar
    Ub53bar = Tb53bar
    Ua1 = -Ud
    Ua2 = Ua
    Ub1 = Ub23
    Ub2 = -Ub53
    U5 = y13*(y1*Uq5 + y3*u3)
    Up5 = y13*(y3*Uq5 - y1*u3)
    Ua1bar = -Udbar
    Ua2bar = Uabar
    Ub1bar = Ub23bar
    Ub2bar = -Ub53bar
    U5bar = y13*(y1*Uq5bar + y3*u3bar)
    Up5bar = y13*(y3*Uq5bar - y1*u3bar)

    Uc = np.transpose([-r*(Ub1+Ua2),1j*r*(Ub1 - Ua2),r*(Ub2 - Ua1),1j*r*(Ub2+Ua1),U5,U6])
    Ucbar = np.array([-r*(Ub1bar+Ua2bar),-1j*r*(Ub1bar - Ua2bar),r*(Ub2bar - Ua1bar),-1j*r*(Ub2bar+Ua1bar),U5bar,U6bar])
    Upc5 = np.transpose([0,0,0,0,Up5,0])
    Upc5bar = np.array([0,0,0,0,Up5bar,0])
    Ub = np.transpose([0,0,0,0,b,0])
    Ubbar = np.array([0,0,0,0,bbar,0])
    UBbar = np.array([0,0,0,0,Bbar,0])
    Ucc = np.transpose([0,0,0,0,c,0])
    Uccbar = np.array([0,0,0,0,cbar,0])
    Uss = np.transpose([0,0,0,0,s,0])
    Ussbar = np.array([0,0,0,0,sbar,0])
    Udd = np.transpose([0,0,0,0,d,0])
    Uddbar = np.array([0,0,0,0,dbar,0])
    Uuu = np.transpose([0,0,0,0,u,0])
    Uuubar = np.array([0,0,0,0,ubar,0])
    UB = np.transpose([0,0,0,0,Bp,0])

    DUc = 1j*gB*a21*(TL1 @ Uc) + 1j*gB*a22*(TL2 @ Uc) + 1j*gB*a23*(TL3 @ Uc) + 1j*gp*b3*((TR3 @ Uc)-(2/3)*(Uno6 @ Uc))
    DUpc5 = 1j*gp*b3*(-2/3)*Upc5
    DUcb = 1j*gp*b3*(1/3)*Ub

    L1i = 1j*(Ucbar @ DUc)
    L2i = 1j*(Upc5bar @ DUpc5)
    L3i = 1j*(Ubbar @ DUcb)
    Lcc = 1j*(Uccbar @ DUc)
    Ldd = 1j*(Uddbar @ DUc)
    Luu = 1j*(Uuubar @ DUc)
    Lss = 1j*(Ussbar @ DUc)
    
    Li = L1i + L2i + L3i + Lcc + Ldd + Luu + Lss

    pl, pr, ga_u = symbols('pl pr ga_u')

    L = 1j*(Li*pl + Ld*pr)*ga_u

    return L

def vector_vertex(lag, qbar, field, q):
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
    fields = ["Z","Zp","WM","Wm","WpM","Wpm","gamma"]
    value_q = [0, 0, 0, 0, 0, 0, 0]
    value_qbar = [0, 0, 0, 0, 0, 0, 0]
    value_field = [0, 0, 0, 0, 0, 0, 0]
    i = 0
    for fila in quark:
        #print(fila,i,end=' ')
        i = i + 1
        if q1 == quark[i-1]:
            indice_q = i-1
            value_q[indice_q] = q
           # print("\n indice: ",i-1,value_q)
        else:
            continue

    j = 0
    for fila in quarkbar:
        j = j + 1
        if qbar1 == quarkbar[j - 1]:
            indice_qbar = j-1
            value_qbar[indice_qbar] = qbar
            #print("\n indice: ", j - 1,value_qbar)
        else:
            continue
    k = 0
    for fila in fields:
        k = k + 1
        if field1 == fields[k - 1]:
            indice_fields = k-1
            value_field[indice_fields] = field
            #print("\n indice: ", k - 1,value_field)
        else:
            continue
    # evaluates the lagrangian
    Lagrangiano_evaluado = lag.subs({d: 0, dbar: 0, u: 0, ubar: 0, s:0, sbar: 0, c: 0, cbar: 0, b: 0,
                                     bbar: 0, t: value_q[0], T: value_q[1], T5: value_q[2], T6: value_q[3],
                                     Tb23: value_q[4], Tb53: value_q[5], Bp: value_q[6], tbar: value_qbar[0],
                                     Tbar: value_qbar[1], T5bar: value_qbar[2], T6bar: value_qbar[3],
                                     Tb23bar: value_qbar[4], Tb53bar: value_qbar[5], Bbar: value_qbar[6],
                                     Z: value_field[0], Zp: value_field[1], WM: value_field[2],
                                     Wm: value_field[3], WpM: value_field[4], Wpm: value_field[5],
                                     gamma: value_field[6]})
    #print(Lagrangiano_evaluado) # se muestra el lagrangiano sólo con: qbar,field,q
    Lag_expan = expand(Lagrangiano_evaluado) # se debe expandir cada paso o lo hace mal
    expr_q = collect(Lag_expan, q).coeff(q, 1)
    expr_q_expan = expand(expr_q)
    expr_qbar = collect(expr_q_expan, qbar).coeff(qbar, 1)
    expr_qbar_expan = expand(expr_qbar)
    expr_field1 = collect(expr_qbar_expan, field).coeff(field, 1)
    expr_field = expand(expr_field1)
    return expr_field


def vector_vertex_value(expr_field, beta_value=1.35, alfa_value=np.pi/2, f_value=1000,
                        y1_value=0.9, y2_value=0.7, y3_value=0.33, cp_value=1, gp_value=0.3528, g_value=0.6528): #, gA: gB):
    '''
    Method to obtain the numerical value of the vertex
    '''
    vertex_value = expr_field.subs({beta: beta_value, alfa: alfa_value, f: f_value, y1: y1_value, y2: y2_value, y3: y3_value, sqrt(2): 1.41421, cp: cp_value, gp: gp_value, g: g_value, gA: gB})
    return vertex_value # return type=sympycore object

    
L_vector = vector_lagrangian()
vertex_expr = vector_vertex(L_vector, "Tbar", "WpM", "Bp")
print(vector_vertex_value(vertex_expr))

