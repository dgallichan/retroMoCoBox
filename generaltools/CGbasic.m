function b = CGbasic(E,m,varargin)
% b = CGbasic(E,m,varargin)
%
% Basic linear conjugate-gradient approximation to inverse

[maxIters, useDisp] = process_options(varargin,'maxIters',15,'useDisp',0);

m = single(m); % data to be reconstructed, [nT, nC], complex
   
cg_a = E'*m;    

outDims = size(cg_a);

b_approx= repmat( single(0+1e-30j), length(cg_a(:)), 1);
p = cg_a(:);
r = cg_a(:);
iteration = 0;

if useDisp
    imab(cg_a)
    title('Conjugate phase recon')
    pause(1)
    tic
end

while iteration < maxIters %&& delta > epsilon %
    
    disp(['Iteration : ' num2str(iteration)])
    
    iteration = iteration + 1;
    
    %%% Linear CG:    
    Ep = E*p;
     
    % Now work out q = E' * Ep;
    q = E' * Ep;

    cg_beta = (r'*r/(p'*q));
    b_approx = b_approx + cg_beta*p; % new estimate of b_approx
    r_new = r - cg_beta*q; % new residual

    p = r_new + (r_new'*r_new/(r'*r))*p; % construct new search direction as linear comb of old direction and new residual

    r = r_new; % replace old r with new r for next iteration    
    
    
    if useDisp
        imab(reshape(b_approx,outDims))
        title(['Iteration: ' num2str(iteration)])
        pause(1)
        toc
        tic
    end
end
    

b = reshape(b_approx,outDims);