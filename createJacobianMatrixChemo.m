function J=createJacobianMatrixChemo(N)

    for i=1:N^2

        % Diffusion terms
        if(mod(i,N)==1), LEFT(i)=i-1+N; else LEFT(i)=i-1; end
        if(mod(i,N)==0), RIGHT(i)=i+1-N; else RIGHT(i)=i+1; end
        if(i>(N-1)*N),   UP(i)=i-(N-1)*N; else UP(i)=i+N; end
        if(i<N+1),       DOWN(i)=(N-1)*N+i; else DOWN(i)=i-N; end

    end
    
    n=1:N^2;            
    a=N^2+1:2*N^2;      
    m=2*N^2+1:3*N^2;    
    c=3*N^2+1:4*N^2;    
    g=4*N^2+1:5*N^2;  
    
    %% Diffusion terms
    n_left=LEFT;
    n_right=RIGHT;
    n_up=UP;
    n_down=DOWN;
    m_left=2*N^2+LEFT;
    m_right=2*N^2+RIGHT;
    m_up=2*N^2+UP;
    m_down=2*N^2+DOWN;
    c_left=3*N^2+LEFT;
    c_right=3*N^2+RIGHT;
    c_up=3*N^2+UP;
    c_down=3*N^2+DOWN;   
    g_left=4*N^2+LEFT;
    g_right=4*N^2+RIGHT;
    g_up=4*N^2+UP;
    g_down=4*N^2+DOWN;    
    
    
    %% n equation  
    n_total_rows=repmat(n,11,1);
    n_total_cols=[n; c; g; n_left; n_right; n_up; n_down; c_left; c_right; c_up; c_down];
    Jn=sparse(n_total_rows, n_total_cols,1,5*N^2,5*N^2);
    
    %% a equation     
    a_total_rows=repmat(a,5,1);
    a_total_cols=[n; a; m; c; g];
    Ja=sparse(a_total_rows, a_total_cols,1,5*N^2,5*N^2);
    
    %% m equation    
    m_total_rows=repmat(m,10,1);
    m_total_cols=[m; c; m_left; m_right; m_up; m_down; c_left; c_right; c_up; c_down];
    Jm=sparse(m_total_rows, m_total_cols,1,5*N^2,5*N^2); 
    
    %% c equation    
    c_total_rows=repmat(c,7,1);
    c_total_cols=[n; a; c; c_left; c_right; c_up; c_down];
    Jc=sparse(c_total_rows, c_total_cols,1,5*N^2,5*N^2);     
    
    
    %% g equation    
    g_total_rows=repmat(g,7,1);
    g_total_cols=[m; a; g; g_left; g_right; g_up; g_down];
    Jg=sparse(g_total_rows, g_total_cols,1,5*N^2,5*N^2);      
    
    %% Add all the bits together
    J=Jn+Ja+Jm+Jc+Jg;
end