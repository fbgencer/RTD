x1,x2 = symbols('x1 x2',real=True)
k1,k2 = symbols('k1 k2')
A1,A2,B1,B2 = symbols('A1 A2 B1 B2')

#M1 = Matrix([[exp(I*k1*x1),exp(-I*k1*x1)],[I*k1*exp(I*k1*x1), -I*k1*exp(-I*k1*x1)]])
M1 = Matrix([[ exp(I*conjugate(k1)*x1), exp(-I*conjugate(k1)*x1) ],[conjugate(-I*k1)*exp(I*conjugate(k1)*x1), -conjugate(-I*k1)*exp(-I*conjugate(k1)*x1) ]])



#M2 = Matrix([[exp(I*k2*x1),exp(-I*k2*x1)],[I*k2*exp(I*k2*x1), -I*k2*exp(-I*k2*x1)]])

M2 = Matrix([[ exp(I*conjugate(k2)*x1),exp(-I*conjugate(k2)*x1) ],[conjugate(-I*k2)*exp(I*conjugate(k2)*x1), -conjugate(-I*k2)*exp(-I*conjugate(k2)*x1) ]])


⎡                 ⎛__   __⎞                      ⎛__   __⎞⎤
⎢ ⎛__   __⎞  ⅈ⋅x₁⋅⎝k₁ - k₂⎠   ⎛  __   __⎞  -ⅈ⋅x₁⋅⎝k₁ + k₂⎠⎥
⎢ ⎝k₁ + k₂⎠⋅ℯ                 ⎝- k₁ + k₂⎠⋅ℯ               ⎥
⎢ ─────────────────────────   ────────────────────────────⎥
⎢              __                           __            ⎥
⎢            2⋅k₂                         2⋅k₂            ⎥
⎢                                                         ⎥
⎢                  ⎛__   __⎞                  ⎛  __   __⎞ ⎥
⎢⎛  __   __⎞  ⅈ⋅x₁⋅⎝k₁ + k₂⎠  ⎛__   __⎞  ⅈ⋅x₁⋅⎝- k₁ + k₂⎠ ⎥
⎢⎝- k₁ + k₂⎠⋅ℯ                ⎝k₁ + k₂⎠⋅ℯ                 ⎥
⎢───────────────────────────  ─────────────────────────── ⎥
⎢              __                           __            ⎥
⎣            2⋅k₂                         2⋅k₂            ⎦
