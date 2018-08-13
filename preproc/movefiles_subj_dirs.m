clear

patients = {'S06',...
            'S07','S08','S09','S10','S11','S12',...
            'S13','S14','S15','S16','S17','S18',...
            'S19','S20','S21','S22','S23','S24'};

%                     'S25','S26','S27','S28','S29','S30',...
%             'S31','S32','S33','S34','S35','S36',... 
%             'S37','S38','S39','S40','S41','S42',...
%             'S43','S44','S45','S46','S47','S48',...
%             'S49','S50','S51','S52','S53','S54'};
        
in_path  = 'D:\ShapeTapper-Analysis\';
out_path = 'E:\ShapeTapper-Analysis\';

for p=1:length(patients)
    
   figs = [in_path patients{p} '\figures\'];
   obs  = [in_path patients{p} '\observed_touchpoints\'];
   
   out_path_patient = [out_path patients{p} '\'];
   
   movefile figs out_path_patient

end %patient loop