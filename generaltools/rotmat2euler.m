function out = rotmat2euler(rotMat)

nT = size(rotMat,3);

out = zeros(3,nT);


% % take the definition from the Wikipedia page:
% % http://en.wikipedia.org/wiki/Rotation_formalisms_in_three_dimensions#Conversion_formulae_between_formalisms
% % Nov 7, 2012
% out(1) = atan2(rotMat(3,1),rotMat(3,2))*180/pi;
% out(2) = acos(rotMat(3,3))*180/pi;
% out(3) = -atan2(rotMat(1,3),rotMat(2,3))*180/pi;


% from /usr/local/fsl/src/miscmaths/miscmaths.cc
% int rotmat2euler(ColumnVector& angles, const Matrix& rotmat)
%     {
%       // uses the convention that R = Rx.Ry.Rz
%       Tracer tr("rotmat2euler");
%       float cz, sz, cy, sy, cx, sx;
%       cy = std::sqrt(Sqr(rotmat(1,1)) + Sqr(rotmat(1,2)));
%       if (cy < 1e-4) {
% 	//cerr << "Cos y is too small - Gimbal lock condition..." << endl;
% 	cx = rotmat(2,2);
% 	sx = -rotmat(3,2);
% 	sy = -rotmat(1,3);
% 	angles(1) = atan2(sx,cx);
% 	angles(2) = atan2(sy,(float)0.0);
% 	angles(3) = 0.0;
%       } else {
% 	// choose by convention that cy > 0
% 	// get the same rotation if: sy stays same & all other values swap sign
% 	cz = rotmat(1,1)/cy;
% 	sz = rotmat(1,2)/cy;
% 	cx = rotmat(3,3)/cy;
% 	sx = rotmat(2,3)/cy;
% 	sy = -rotmat(1,3);
% 	//atan2(sin,cos)  (defined as atan2(y,x))
% 	angles(1) = atan2(sx,cx);
% 	angles(2) = atan2(sy,cy);
% 	angles(3) = atan2(sz,cz);
%       }
%       return 0;
%     }
%     


for iT = 1:nT
    
    cy = sqrt(rotMat(1,1,iT).^2 + rotMat(1,2,iT).^2);
    if cy < 1e-4
        cx = rotMat(2,2,iT);
        sx = -rotMat(3,2,iT);
        sy = -rotMat(1,3,iT);
        out(1,iT) = atan2(sx,cx);
        out(2,iT) = atan2(sy,0);
        out(3,iT) = 0;
    else
        cz = rotMat(1,1,iT)/cy;
        sz = rotMat(1,2,iT)/cy;
        cx = rotMat(3,3,iT)/cy;
        sx = rotMat(2,3,iT)/cy;
        sy = -rotMat(1,3,iT);
        out(1,iT) = atan2(sx,cx);
        out(2,iT) = atan2(sy,cy);
        out(3,iT) = atan2(sz,cz);
    end
    
end

out = out * 180/pi; 
    