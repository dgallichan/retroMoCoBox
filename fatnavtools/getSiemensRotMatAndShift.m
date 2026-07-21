function out = getSiemensRotMatAndShift(hdr)
% input 'hdr' should be twix_obj.hdr
%
% matrix O corresponds to the Siemens internal usage, such that( PE ; RO ; Slice) = O *(x_anat - s_anat)
% Anat is Sag Cor Trans = Left Posterior Head
%
% Assumption for checks : Siemens uses FFT for kspace -> image space ('confirmed' by IDEA Discussion board)
%
% Adapted from code by Frederic Gretsch (gretsch@epfl.ch) by gallichand@cardiff.ac.uk 
% output in ColLinPar = RO PE SLICE 
%

out = struct();

sNormal = hdr.MeasYaps.sSliceArray.asSlice{1}.sNormal;
if ~isfield(sNormal,'dCor'),  sNormal.dCor = 0; end
if ~isfield(sNormal,'dTra'),  sNormal.dTra = 0; end
if ~isfield(sNormal,'dSag'),  sNormal.dSag = 0; end

if ~isfield(hdr.MeasYaps.sSliceArray.asSlice{1},'sPosition')
    out.sPosition.dCor = 0;
    out.sPosition.dTra = 0;
    out.sPosition.dSag = 0;
else
    out.sPosition = hdr.MeasYaps.sSliceArray.asSlice{1}.sPosition;
    if ~isfield(out.sPosition,'dCor'), out.sPosition.dCor = 0; end
    if ~isfield(out.sPosition,'dTra'), out.sPosition.dTra = 0; end
    if ~isfield(out.sPosition,'dSag'), out.sPosition.dSag = 0; end
end


% Sag x Cor = Tra
nSag= sNormal.dSag ;
nCor= sNormal.dCor ;
nTra= sNormal.dTra ;

n=[nSag nCor nTra].';

if isfield(hdr.MeasYaps.sSliceArray.asSlice{1},'dInPlaneRot')
    dInPlaneRot = hdr.MeasYaps.sSliceArray.asSlice{1}.dInPlaneRot;
else
    dInPlaneRot = 0;
end

SAGITTAL=0;
CORONAL=1;
TRANSVERSE=2;
[dGp ,dGr,lCase ] = fGSLCalcPRS (n,dInPlaneRot    );

switch lCase
    case 0
        out.baseOrientation = 'Sagittal';
    case 1
        out.baseOrientation = 'Coronal';
    case 2
        out.baseOrientation = 'Transverse';
end
    

dGr=dGr(:);
dGp=dGp(:);
%         Now we have the three vectors in patient coords

% Phase RO Slice to anat
U = [dGp,dGr,n];

O=U';

% in Col Lin Par 
out.RotMat = [0 1 0; 1 0 0 ; 0 0 1 ] * O;
out.dInPlaneRot = dInPlaneRot;
out.sNormal = sNormal;
out.Shifts_SagCorTra = [out.sPosition.dSag out.sPosition.dCor out.sPosition.dTra];


function b= fGSLAlmEqual ...
            ( ...
            dArg1,                ...         /* IMP: first argument              */
            dArg2                 ...         /* IMP: second argument             */
            )
        dTmp = dArg1 - dArg2;
        b= (dTmp >= -1.e-6  &&  dTmp <= 1.e-6);
end




function lCase = fGSLClassOri( ...
            dSagComp  ,...   % normal vector components
            dCorComp , ...
            dTraComp ...
            )
        
        %   //---------------------------------------------------------------------------
        %   // For performance reasons calculate some tmp values
        %   //---------------------------------------------------------------------------
        dAbsSagComp     = abs(dSagComp);
        dAbsCorComp     = abs(dCorComp);
        dAbsTraComp     = abs(dTraComp);
        bAlmEqualSagCor = fGSLAlmEqual(dAbsSagComp, dAbsCorComp);
        bAlmEqualSagTra = fGSLAlmEqual(dAbsSagComp, dAbsTraComp);
        bAlmEqualCorTra = fGSLAlmEqual(dAbsCorComp, dAbsTraComp);
        
        %   //---------------------------------------------------------------------------
        %   // Check all values to determine the slice orientation (sag, cor, tra)
        %   //---------------------------------------------------------------------------
        if ((bAlmEqualSagCor              &&  bAlmEqualSagTra)             || ...
                (bAlmEqualSagCor              &&  (dAbsSagComp < dAbsTraComp)) || ...
                (bAlmEqualSagTra              &&  (dAbsSagComp > dAbsCorComp)) || ...
                (bAlmEqualCorTra              &&  (dAbsCorComp > dAbsSagComp)) || ...
                ((dAbsSagComp > dAbsCorComp)  &&  (dAbsSagComp < dAbsTraComp)) || ...
                ((dAbsSagComp < dAbsCorComp)  &&  (dAbsCorComp < dAbsTraComp)) || ...
                ((dAbsSagComp < dAbsTraComp)  &&  (dAbsTraComp > dAbsCorComp)) || ...
                ((dAbsCorComp < dAbsTraComp)  &&  (dAbsTraComp > dAbsSagComp)))
            
            %     //-------------------------------------------------------------------------
            %     // Mainly transverse...
            %     //-------------------------------------------------------------------------
            lCase = TRANSVERSE;
            
        elseif ((bAlmEqualSagCor              &&  (dAbsSagComp > dAbsTraComp)) || ...
                (bAlmEqualSagTra              &&  (dAbsSagComp < dAbsCorComp)) || ...
                ((dAbsSagComp < dAbsCorComp)  &&  (dAbsCorComp > dAbsTraComp)) || ...
                ((dAbsSagComp > dAbsTraComp)  &&  (dAbsSagComp < dAbsCorComp)) || ...
                ((dAbsSagComp < dAbsTraComp)  &&  (dAbsTraComp < dAbsCorComp)))
            
            %     //-------------------------------------------------------------------------
            %     // Mainly coronal...
            %     //-------------------------------------------------------------------------
            lCase = CORONAL;
            
        elseif ((bAlmEqualCorTra              &&  (dAbsCorComp < dAbsSagComp)) || ...
                ((dAbsSagComp > dAbsCorComp)  &&  (dAbsSagComp > dAbsTraComp)) || ...
                ((dAbsCorComp > dAbsTraComp)  &&  (dAbsCorComp < dAbsSagComp)) || ...
                ((dAbsCorComp < dAbsTraComp)  &&  (dAbsTraComp < dAbsSagComp)))
            
            %     //-------------------------------------------------------------------------
            %     // Mainly sagittal...
            %     //-------------------------------------------------------------------------
            lCase = SAGITTAL;
            
        else
            
            %     //-------------------------------------------------------------------------
            %     // Invalid slice orientation...
            %     //-------------------------------------------------------------------------
            error('HUH')
        end
        
    end

    function [dGp ,dGr,lCase ] = fGSLCalcPRS ...
            ( ...
            dGs,     ...       
            dPhi        ...     
            )
        
        %/* Orientation (SAGITTAL, CORONAL or TRANSVERSE)  */
        
        lCase = fGSLClassOri (dGs(SAGITTAL+1), dGs(1+CORONAL), dGs(1+TRANSVERSE));
        
        
        switch (lCase)
            
            case TRANSVERSE
                dGp(1) = 0.;
                dGp(2) = dGs(3) * sqrt (1. / (dGs(2) * dGs(2) + dGs(3) * dGs(3)));
                dGp(3) = -dGs(2) * sqrt (1. / (dGs(2) * dGs(2) + dGs(3) * dGs(3)));
                
            case CORONAL
                dGp(1) = dGs(2) * sqrt (1. / (dGs(1) * dGs(1) + dGs(2) * dGs(2)));
                dGp(2) = -dGs(1) * sqrt (1. / (dGs(1) * dGs(1) + dGs(2) * dGs(2)));
                dGp(3) = 0.;
                
            case SAGITTAL
                dGp(1) = -dGs(2) * sqrt (1. / (dGs(1) * dGs(1) + dGs(2) * dGs(2)));
                dGp(2) = dGs(1) * sqrt (1. / (dGs(1) * dGs(1) + dGs(2) * dGs(2)));
                dGp(3) = 0.;
                
        end
        
        %   /*------------------------------------------------------------------------*/
        %   /*  Calculate GR = GS x GP                                                */
        %   /*------------------------------------------------------------------------*/
        dGr(1) = dGs(2) * dGp(3) - dGs(3) * dGp(2);
        dGr(2) = dGs(3) * dGp(1) - dGs(1) * dGp(3);
        dGr(3) = dGs(1) * dGp(2) - dGs(2) * dGp(1);
        
        
        if (dPhi ~= 0)
            
            %     /*----------------------------------------------------------------------*/
            %     /* Rotate around the S axis                                             */
            %     /*----------------------------------------------------------------------*/
            
            dGp(1) = cos (dPhi) * dGp(1) - sin (dPhi) * dGr(1);
            dGp(2) = cos (dPhi) * dGp(2) - sin (dPhi) * dGr(2);
            dGp(3) = cos (dPhi) * dGp(3) - sin (dPhi) * dGr(3);
            
            %     /*----------------------------------------------------------------------*/
            %     /* Calculate new GR = GS x GP                                           */
            %     /*----------------------------------------------------------------------*/
            dGr(1) = dGs(2) * dGp(3) - dGs(3) * dGp(2);
            dGr(2) = dGs(3) * dGp(1) - dGs(1) * dGp(3);
            dGr(3) = dGs(1) * dGp(2) - dGs(2) * dGp(1);
            
        end
        
    end

    end

