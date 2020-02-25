fx = jacobian(FF,XX);
fxp = jacobian(FF,XXp);
fy = jacobian(FF,YY);
fyp = jacobian(FF,YYp);
fz = jacobian(FF,ZZ);
fzp = jacobian(FF,ZZp);

fzz = reshape(jacobian(fz(:),ZZ),nx+ny,nz,nz);
fxz = reshape(jacobian(fx(:),ZZ),nx+ny,nx,nz);
fxx = reshape(jacobian(fx(:),XX),nx+ny,nx,nx);
fyz = reshape(jacobian(fy(:),ZZ),nx+ny,ny,nz);
fyx = reshape(jacobian(fy(:),XX),nx+ny,ny,nx);
fyy = reshape(jacobian(fy(:),YY),nx+ny,ny,ny);

fzpzp = reshape(jacobian(fzp(:),ZZp),nx+ny,nz,nz);
fxpzp = reshape(jacobian(fxp(:),ZZp),nx+ny,nx,nz);
fxpxp = reshape(jacobian(fxp(:),XXp),nx+ny,nx,nx);
fypzp = reshape(jacobian(fyp(:),ZZp),nx+ny,ny,nz);
fypxp = reshape(jacobian(fyp(:),XXp),nx+ny,ny,nx);
fypyp = reshape(jacobian(fyp(:),YYp),nx+ny,ny,ny);