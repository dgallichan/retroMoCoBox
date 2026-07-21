function [b, b_iters] = CGbasic(E,m,varargin)
% b = CGbasic(E,m,varargin)
%
% Basic linear conjugate-gradient approximation to inverse
%
% updated Sep 25 to attempt to add 'hot start' (not yet correct!)

[maxIters, useDisp, b_start] = process_options(varargin,'maxIters',15,'useDisp',0,'b_start',[]);

m = single(m); % data to be reconstructed, [nT, nC], complex
   
cg_a = E'*m;    

outDims = size(cg_a);

if nargout > 1
    b_iters = zeros([prod(outDims) maxIters]);
end

if isempty(b_start)
    b_approx= repmat( single(0+1e-30j), length(cg_a(:)), 1);
    p = cg_a(:);
    r = cg_a(:);
else
    % this bit is probably not right yet!
    b_approx = b_start(:);
    p = b_start(:);
    r = cg_a(:);
end

iteration = 0;

if useDisp
    imab(cg_a)
    title('Conjugate phase recon')
    pause(1)
    tic
end

while iteration < maxIters %&& delta > epsilon %
    
    % disp(['Iteration : ' num2str(iteration)])
    
    iteration = iteration + 1;
    
    %%% Linear CG:    
    Ep = E*p;
     
    % Now work out q = E' * Ep;
    q = E' * Ep;

    cg_beta = (r'*r/(p'*q(:)));
    b_approx = b_approx + cg_beta*p; % new estimate of b_approx

    if nargout > 1
        b_iters(:,iteration) = b_approx;
    end

    r_new = r - cg_beta*q(:); % new residual

    p = r_new + (r_new'*r_new/(r'*r))*p; % construct new search direction as linear comb of old direction and new residual

    r = r_new; % replace old r with new r for next iteration    
    
    
    if useDisp == 2
        imab(reshape(b_approx,outDims))
        title(['Iteration: ' num2str(iteration)])
        pause(1)
    elseif useDisp == 1
        toc
        tic
    end
end
    

b = reshape(b_approx,outDims);
if nargout > 1
    b_iters = reshape(b_iters,[outDims maxIters]);
end