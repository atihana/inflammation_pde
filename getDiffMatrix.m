function DMAT=getDiffMatrix(n,h,bcType)

DMAT=zeros(n^2);

switch bcType
    
    case 'Periodic'
        
        for i=1:n^2
            
            %Periodic BCs
            if(mod(i,n)==1), LEFT=i-1+n; else LEFT=i-1; end
            if(mod(i,n)==0), RIGHT=i+1-n; else RIGHT=i+1; end
            if(i>(n-1)*n),   UP=i-(n-1)*n; else UP=i+n; end
            if(i<n+1),       DOWN=(n-1)*n+i; else DOWN=i-n; end
            
            DMAT(i,i)=-4;
            DMAT(i,UP)=1;
            DMAT(i,DOWN)=1;
            DMAT(i,LEFT)=1;
            DMAT(i,RIGHT)=1;
            
        end
        
    case 'Dirichlet'
        
        for i=1:n^2
            
            %Zero Dirichlet BC
            
            if(~(mod(i,n)==1)), DMAT(i,i-1)=1; end
            if(~(mod(i,n)==0)), DMAT(i,i+1)=1; end
            if(~(i>(n-1)*n)),   DMAT(i,i+n)=1; end
            if(~(i<n+1)),       DMAT(i,i-n)=1; end
            if(mod(i,n)==1 || mod(i,n)==0 || i>(n-1)*n || i<n+1)
                DMAT(i,:)=0;
                %DMAT(i,i)=1;
            else
                DMAT(i,i)=-4;
            end
            
        end
        
    case 'Neumann'
        
        for i=1:n^2
            
            %Neumann BCs
            if(mod(i,n)==1)
                DMAT(i,i+1)=2;
            elseif(mod(i,n)==0)
                DMAT(i,i-1)=2;
            else
                DMAT(i,i+1)=1;
                DMAT(i,i-1)=1;
            end
            if(i>(n-1)*n)
                DMAT(i,i-n)=2;
            elseif(i<n+1)
                DMAT(i,i+n)=2;
            else
                DMAT(i,i+n)=1;
                DMAT(i,i-n)=1;
            end
            DMAT(i,i)=-4;         
            
        end
        
        %DMAT=DMAT/h^2;
        
end