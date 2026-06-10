# === Validation against supervisor's Maple results ===
print("\n=== VALIDATION: n=10, beta=0.09, c=0.99 ===")
n_val, beta_val, c_val = 10, 0.09, 0.99

zeros = meixner_zeros(n_val, beta_val, c_val)
x1n = zeros[0]
xnn = zeros[-1]
print(f"x_{{1,10}} = {x1n:.15f}  (Maple: 0.886751385278167)")
print(f"x_{{10,10}} = {xnn:.15f}  (Maple: 2814.03598753338)")

# k=0: F_n^+(0) should match Maple's boundMeix1 = 2555.231...
vals0 = Fn0(n_val, beta_val, c_val)
Fplus0  = max(vals0)
Fminus0 = min(vals0)
print(f"\nk=0:")
print(f"F_n^+(0)  = {Fplus0:.15f}  (Maple: 2555.23118013068)")
print(f"F_n^-(0)  = {Fminus0:.15f}")

# k=6: F_n^-(6) should match Maple's boundMeix2 = 0.8867...
vals6 = Fn6(n_val, beta_val, c_val)
Fplus6  = max(vals6)
Fminus6 = min(vals6)
print(f"\nk=6:")
print(f"F_n^+(6)  = {Fplus6:.15f}")
print(f"F_n^-(6)  = {Fminus6:.15f}  (Maple: 0.886768650875350)")

# Check inequalities
print(f"\nCheck x_{{1,n}} < F_n^-(0): {x1n:.6f} < {Fminus0:.6f} => {x1n < Fminus0}")
print(f"Check F_n^+(0) < x_{{n,n}}: {Fplus0:.6f} < {xnn:.6f} => {Fplus0 < xnn}")
print(f"Check x_{{1,n}} < F_n^-(6): {x1n:.6f} < {Fminus6:.6f} => {x1n < Fminus6}")
print(f"Check F_n^+(6) < x_{{n,n}}: {Fplus6:.6f} < {xnn:.6f} => {Fplus6 < xnn}")
