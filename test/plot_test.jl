using Ellipsoids

c = [1.5; 1.5]
P = [4.0 0.5;       
     0.5 6.0]
a = 2.5 #1.52
c0 = [1.6+a; 1.4+a]
P0 = [0.4 -0.1;
     -0.1 0.5]

El = Ellipsoids.Ellipsoid(P, c)
El0 = Ellipsoids.Ellipsoid(P0, c0)
println(Ellipsoids.intersect(El,El0))

ell_ast, β_ast =  Ellipsoids.get_ℓ_ast_intersect(El, El0)
Elnew = El0*(ell_ast)#Ellipsoids.Ellipsoid(P0, ell_ast*c0)

#plot()
Ellipsoids.plot_ellipsoid!(El, label="El")
Ellipsoids.plot_ellipsoid!(El0, label="El0")
Ellipsoids.plot_ellipsoid!(Elnew, label="Elnew")
