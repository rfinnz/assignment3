%% Part 2 
% in this part of the assignment, the code I completed in assignment 2 was
% used to calculate the potential with the bottle-neck inserted. From this
% potential, the electric field was calculated using quiver.

%setting up dimensions and matrices
nx = 50;                
ny = (3/2)*nx;            
G = sparse(nx*ny);      
Op = zeros(1, nx*ny);    


Sigmatrix = zeros(nx, ny);    % a sigma matrix is required for this part
Sig1 = 1;                     % sigma value given outside the box
Sig2 = 10^-2;                 % sigma value given inside the box

%The box will be difined using a 1x4 matrix containing it's dimensions
box = [nx*2/5 nx*3/5 ny*2/5 ny*3/5]; 

for i = 1:nx
    
    for j = 1:ny
        
        
        if i > box(1) && i < box(2) && (j < box(3)||j > box(4))
            Sigmatrix(i, j) = Sig2;
            
        else
            Sigmatrix(i, j) = Sig1;
            
            
        end
    end
end

% Filling in G matrix with corresponding bottleneck conditions
for x = 1:nx
    
    for y = 1:ny
        
        n = y + (x-1)*ny;
        nposx = y + (x+1-1)*ny;
        nnegx = y + (x-1-1)*ny;
        nposy = y + 1 + (x-1)*ny;
        nnegy = y - 1 + (x-1)*ny;
        
        if x == 1
            
            G(n, :) = 0;
            G(n, n) = 1;
            Op(n) = 1;
            
        elseif x == nx
            
            G(n, :) = 0;
            G(n, n) = 1;
            Op(n) = 0;
            
        elseif y == 1

            G(n, nposx) = (Sigmatrix(x+1, y) + Sigmatrix(x,y))/2;
            G(n, nnegx) = (Sigmatrix(x-1, y) + Sigmatrix(x,y))/2;
            G(n, nposy) = (Sigmatrix(x, y+1) + Sigmatrix(x,y))/2;            
            G(n, n) = -(G(n,nposx)+G(n,nnegx)+G(n,nposy));
            
        elseif y == ny
            
            G(n, nposx) = (Sigmatrix(x+1, y) + Sigmatrix(x,y))/2;
            G(n, nnegx) = (Sigmatrix(x-1, y) + Sigmatrix(x,y))/2;
            G(n, nnegy) = (Sigmatrix(x, y-1) + Sigmatrix(x,y))/2;
            G(n, n) = -(G(n,nposx)+G(n,nnegx)+G(n,nnegy));
            
        else
            
            G(n, nposx) = (Sigmatrix(x+1, y) + Sigmatrix(x,y))/2;
            G(n, nnegx) = (Sigmatrix(x-1, y) + Sigmatrix(x,y))/2;
            G(n, nposy) = (Sigmatrix(x, y+1) + Sigmatrix(x,y))/2;
            G(n, nnegy) = (Sigmatrix(x, y-1) + Sigmatrix(x,y))/2;
            G(n, n) = -(G(n,nposx)+G(n,nnegx)+G(n,nposy)+G(n,nnegy));
            
        end
    end
end


%Voltage matrix calculation
Voltage = G\Op';


sol = zeros(ny, nx, 1);
for i = 1:nx
    for j = 1:ny
        n = j + (i-1)*ny;
        sol(j,i) = Voltage(n);
    end
end

%V(x,y) Surface Plot
figure(1)
surf(sol)
axis tight
xlabel("X position")
ylabel("Y position")
zlabel("Voltage")
view([40 30]);
title("Voltage Surface Plot with Given Bottleneck Conditions")

%The electric field can be derived from the surface voltage using a
%gradient

[elecx, elecy] = gradient(sol);

%plotting the electric field from the potential using quiver

figure(2)
quiver(-elecx, -elecy);
axis tight
title("2-D electric field vector plot")

