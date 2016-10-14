%Lennard Jones Water model for TIP3P

function [ff,u,config] = fcalc_2(xx,boxl,np)
  % function [ff,u] = fcalc(xx,nn,ns,p,rc,boxl)
  %   Computation of interfarticle forces using
  %   Lennard-Jones interatomic potential
  %   xx --> nx3 array of particle coordinates
  %   boxl --> size of (cubic) simulation box
  %   np --> number of particles
epsilonHH=0.0460;
epsilonOO=0.1521;
sigmaHH=.4;
sigmaOO=.0460;
fc=zeros(np,3);
ff = zeros(np,3);
p = zeros(np,1);
dr= zeros(np,3);
config=0;

mol=size(xx,1)/3;
    
  for i=1:np-1
      
      for j=i+1:np 
                  
        dr(i,:)=xx(i,:)-xx(j,:);  
        dr(i,:)=dr(i,:)-boxl*round((dr(i,:)/boxl));
        r=sqrt(sum(dr(i,:).^2));

        if ( (i<=mol && j<=mol) || (i>=mol*2+1 && j>=mol*2+1) ||(i<=mol && j>=mol*2+1))
              epsilon=epsilonHH;
              sigma=sigmaHH;
        elseif ( (i<=mol && j>=mol+1 && j<=2*mol) || (i>=mol+1 && i<=mol*2 && j>=2*mol+1))
              epsilon=sqrt(epsilonHH*epsilonOO);
              sigma=0.5*(sigmaHH+sigmaOO);
        elseif (j>=mol+1 && j<=2*mol && i>=mol+1 && i<=mol*2)
              epsilon=epsilonOO;
              sigma=sigmaOO;
        else
            'error: epsilon sigma!';
        end
       
          if(r>2.5*sigma)    
          p(i) =p(i)+4*epsilon*((sigma/r)^12 - (sigma/r)^6);
          %u(j)=u(j)-u(i);
          f=4*epsilon*(-12*sigma^12*r^(-13)+6*sigma^6*r^-7);
          fc(i,:)=f*dr(i,:)./r;
          ff(i,:)=ff(i,:)+fc(i,:);
          ff(j,:)=ff(j,:)-fc(i,:);
          config=config+abs(dot(ff(i,:)',dr(i,:)'));
          end
      end
  end
    u=sum(p);
  end
      
      