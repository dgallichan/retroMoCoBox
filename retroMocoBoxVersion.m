function out = retroMocoBoxVersion
out = '1.0.1dev';
% 1.0.1dev - Update how the circshift is handled
% 1.0.0dev - MAJOR OVERHAUL - 
%          Highlights
%             - Orientation of host should use Siemens rules for getting
%               the full coordinate transform and (hopefully!) handle more 
%               orientation possibilities gracefully  
%             - MoCo now performed in native RPS coordinate system of host 
%               acquisition (easier to converse with other recon pipelines)                
%             - BUG FIX - If your isocentre is offset significantly from the
%               centre of the head (seems to be the case for UHF) then the
%               motion parameters are (hopefully!) now calculated
%               appropriately (they definitely were not before!)
%             - Added a 'standard' FatNav (average of 5 subjects) that can
%               be used to define a mask in standard space for where the
%               FatNav can best be assumed to be rigid prior to
%               registration (see Marchetto et al, MRM 2023)
%             - As the FatNavs are often wrapped in the PE direction (nose 
%               looking at back of head, added a feature to detect the
%               minimum signal in PE direction, then circularly shift the
%               FatNav to centre it - and (hopefully!) account for this
%               shift properly in the motion parameter estimation.
% 0.10.1dev - tried to clean-up handling of figure windows which can steal focus
% 0.10.0dev - getting multi-echo GRE to work with ASPIRE coil combination, June 2022
% 0.9.0dev - last major change was fixing euler2rmat.m, October 2019


