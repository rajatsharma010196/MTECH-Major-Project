function N=get_N(p_type,xi)

    N_lib={[1];
    [(1-xi)/2,(1+xi)/2];
    [-xi*(1-xi)/2,(1+xi)*(1-xi),xi*(1+xi)/2];
    };
    
    N_g=N_lib(p_type+1);
    N=cell2mat(N_g);
end
