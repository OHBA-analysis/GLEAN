function a = riemann_dist(A,B)

a = sqrt(sum(log(eig(A,B)).^2));