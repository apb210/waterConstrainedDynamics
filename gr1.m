function [rc,gr] = gr1(xx,nbins,rcut,boxl)
  n = size(xx,1);
  mol=n/3;
  % rcut=10;
  
  % Bin edges
  rr = (0:nbins)/nbins * rcut;

  % Variable for storing histogram of distances
  gr = zeros(1,nbins+1);
  xx=xx(mol:2*mol,:);  
  % Compute all possible distances between particles
  for i=1:mol
    X = [xx(:,1)-xx(i,1) , xx(:,2)-xx(i,2) , xx(:,3)-xx(i,3)];
    X = X - boxl*round(X*(1/boxl));
    X = sqrt(sum(X'.^2));

    % Add distances to particle i to histogram
    gr = gr + histc(X,rr);
  end

  % Return center of bins for presentation
  rc = 0.5*(rr(1:nbins)+rr(2:nbins+1));

  % Comparison with the Ideal gas 
  gr = gr(1:nbins) ./ (4/3*pi*(rr(2:nbins+1).^3-rr(1:nbins).^3))/n;

  % Take away the delta function at r=0
  gr(1) = 0;
