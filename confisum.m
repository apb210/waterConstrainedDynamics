function config2 = confisum(L,OHbl,HHbl,mol)
    config2=0;
    for i=1:mol^3
        config2=config2+(L(1)*OHbl-L(3)*HHbl)*OHbl;
        config2=config2+(L(2)*OHbl-L(1)*OHbl)*OHbl;
        config2=config2+(L(3)*HHbl-L(2)*OHbl)*HHbl;
      
    end
    


end



