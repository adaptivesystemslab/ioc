clear all;clc;close all;

modelName='SQUAT_7DOF';%'STS_frame_04'


fileName = [modelName ];
fId = fopen(['../Symoro/' fileName,'.jac'],'r+')


copyfile(['../Symoro/' fileName,'.jac'],['../Symoro/' fileName,'.m']);

fIdstr = fscanf(fId,'%c');


fIdstr = strrep(fIdstr,'Sin','sin');
fIdstr = strrep(fIdstr,'Cos','cos');
fIdstr = strrep(fIdstr,[' + ' char(13)],[' + ... ' char(13)]);
fIdstr = strrep(fIdstr,[' - ' char(13)],[' - ... ' char(13)]);
fIdstr = strrep(fIdstr,['*' char(13)],['* ... ' char(13)]);
 


global_ind=strfind(fIdstr,'global');
size_ind=strfind(fIdstr,'Size=(');
close_par_ind=strfind(fIdstr,')');
equal_ind=strfind(fIdstr,'=');
J_ind=strfind(fIdstr,'J(');
char13_ind=strfind(fIdstr,char(13));

dcl_out_ind=strfind(fIdstr,'Jacobian Matrix')
eval(['Jnan=nan(' fIdstr(size_ind+6:dcl_out_ind-5) ');' ]);



if size(J_ind)~=size(Jnan,1)*size(Jnan,2)
    disp('Eror jacobian size')
end


NbSegment=11;
Nbq=7;

%delare joint angles and segment length
for ii=1:NbSegment
    eval(['syms L',num2str(ii)]);
end
for ii=1:Nbq
    eval(['syms q',num2str(ii)]);
end


syms J

for i=1:size(Jnan,1)*size(Jnan,2)-1;

%equal_next=equal_ind(find(equal_ind>J_ind(i),1));

eval([fIdstr(J_ind(i):J_ind(i+1)-1) ';']);
%strfind(fIdstr(J_ind(i):J_ind(i+1)-1),[' + ' char(13)])+J_ind(i)
%strfind(fIdstr(J_ind(i):J_ind(i+1)-1),[' - ' char(13)])+J_ind(i)

end


%%
%eval([fIdstr(J_ind(end):char13_ind(find(char13_ind>J_ind(end),1)) ) ';'])


%%
matlabFunction(J,'file',[modelName '_Ext_wrenches_Jacobian.m']);%,'Optimize',false)
fclose(fId)

%read the string form matlab function then delete the file to create a new
%one with new header
fId = fopen([modelName '_Ext_wrenches_Jacobian.m'],'r+')
fIdstr = fscanf(fId,'%c');
fclose(fId)
delete([modelName '_Ext_wrenches_Jacobian.m'])

fId = fopen([modelName '_Ext_wrenches_Jacobian.m'],'w+')


fIdstr = strrep(fIdstr,'function','%function');
new_L_str=[];
for i=1:NbSegment
new_L_str=[new_L_str 'L' num2str(i) '=param.L(' num2str(i) ');' char(13)];
end

new_q_str=[];
for i=1:Nbq
new_q_str=[new_q_str 'q' num2str(i) '=q(' num2str(i) ',:);' char(13)];
end


new_head_str=[];
new_head_str=[new_head_str 'function J=' modelName '_Ext_wrenches_Jacobian(q,param)' char(13) new_L_str char(13) new_q_str char(13) fIdstr]


fId = fopen([modelName '_Ext_wrenches_Jacobian.m'],'w+')
fprintf(fId, '%c',[new_head_str])
fclose(fId);




% 
% %%
% fcte_ind=strfind(fIdstr,'% Function description:');
% 
% space_ind=strfind(fIdstr,' ');
% return_ind=strfind(fIdstr,char(10));
% 
% 
% 
% %Search all the q in input declaration header
% %then creat string for q dq qnd ddq allocation
% q_ind_first=strfind(fIdstr,'q');
% q_ind_first=q_ind_first(find(q_ind_first<dcl_out_ind ));
% q_ind_first=q_ind_first(find(q_ind_first>dcl_inp_ind ));
% 
% new_q_str=[];
% new_dq_str=[];
% new_ddq_str=[];
% 
% 
% j=1;
% for i=1:length(q_ind_first)
%     
%     space_next=space_ind(find(space_ind>q_ind_first(i),1));
%     return_next=return_ind(find(return_ind>q_ind_first(i),1));
%     
%     if fIdstr(q_ind_first(i)-1)=='d'
%         if fIdstr(q_ind_first(i)-2)=='d'
%             if return_next<space_next
%                 new_ddq_str=[new_ddq_str fIdstr(q_ind_first(i)-2:return_next-2) '=ddq(',fIdstr(q_ind_first(i)+1:return_next-2),',:);' char(10)];
%             else
%                 new_ddq_str=[new_ddq_str fIdstr(q_ind_first(i)-2:space_next-1) '=ddq(',fIdstr(q_ind_first(i)+1:space_next-1),',:);' char(10)];
%             end
%         else
%             if return_next<space_next
%                 new_dq_str=[new_dq_str fIdstr(q_ind_first(i)-1:return_next-2) '=dq(',fIdstr(q_ind_first(i)+1:return_next-2),',:);' char(10)];
%             else
%                 new_dq_str=[new_dq_str fIdstr(q_ind_first(i)-1:space_next-1) '=dq(',fIdstr(q_ind_first(i)+1:space_next-1),',:);' char(10)];
%             end
%         end
%         
%     else
%         if return_next<space_next
%             new_q_str=[new_q_str fIdstr(q_ind_first(i):return_next-2) '=q(',fIdstr(q_ind_first(i)+1:return_next-2),',:);' char(10)];
%         else
%             new_q_str=[new_q_str fIdstr(q_ind_first(i):space_next-1) '=q(',fIdstr(q_ind_first(i)+1:space_next-1),',:);' char(10)];
%         end
%     end
%     
% end
% 
% 
% 
% % j      W0     WP0    V0     VP0    G
% 
% % 1      WX0    WPX0   VX0    VPX0   0
% 
% % 2      WY0    WPY0   VY0    VPY0   0
% 
% % 3      WZ0    WPZ0   VZ0    VPZ0   G3
% 
% %%Create string for base-link velocity and accelerations and gravity
% 
% 
% new_base0_str=['G2=-9.81;WX0=base0.WX0;WY0=base0.WY0;WZ0=base0.WZ0; WPX0=base0.WPX0;WY0=base0.WPY0;WPZ0=base0.WPZ0; VX0=base0.VX0;VY0=base0.VY0;VZ0=base0.VZ0; VPX0=base0.VPX0;VPY0=base0.VPY0;WZ0=base0.VPZ0; '];
% 
% 
% 
% %Search all the MX in input declaration header then allocate the inertial
% %parameters.
% MX_ind_first=strfind(fIdstr,'MX');
% MX_ind_first=MX_ind_first(find(MX_ind_first<dcl_out_ind ));
% MX_ind_first=MX_ind_first(find(MX_ind_first>dcl_inp_ind ));
% 
% new_COM_str=[];
% new_Mass_str=[];
% new_Inertia_str=[];
% new_FS_str=[];
% new_FV_str=[];
% new_F_str=[];
% new_C_str=[];
% 
% for i=1:length(MX_ind_first)
%     new_Mass_str=[new_Mass_str 'M' num2str(i) '=param.M(' num2str(i) ');'  char(10)];
%     new_COM_str=[new_COM_str 'MX' num2str(i) '=param.MX(' num2str(i) ');' 'MY' num2str(i) '=param.MY(' num2str(i) ');' 'MZ' num2str(i) '=param.MZ(' num2str(i) ');' char(10)];
%     new_Inertia_str=[new_Inertia_str 'XX' num2str(i) '=param.XX(' num2str(i) ');' 'XY' num2str(i) '=param.XY(' num2str(i) ');' 'XZ' num2str(i) '=param.XZ(' num2str(i) ');' char(10)];
%     new_inertia_str=[new_Inertia_str 'YX' num2str(i) '=param.YX(' num2str(i) ');' 'YY' num2str(i) '=param.YY(' num2str(i) ');' 'YZ' num2str(i) '=param.YZ(' num2str(i) ');' char(10)];
%     new_Inertia_str=[new_Inertia_str 'ZX' num2str(i) '=param.ZX(' num2str(i) ');' 'ZY' num2str(i) '=param.ZY(' num2str(i) ');' 'ZZ' num2str(i) '=param.ZZ(' num2str(i) ');' char(10)];
%     new_FS_str=[new_FS_str 'FS' num2str(i) '=param.FS(' num2str(i) ',:);'  char(10)];
%     new_FV_str=[new_FS_str 'FV' num2str(i) '=param.FV(' num2str(i) ',:);'  char(10)];
%     new_F_str=[new_F_str 'FX' num2str(i) '=ext_wrenches.FX(' num2str(i) ',:);' 'FY' num2str(i) '=param.FY(' num2str(i) ',:);' 'FZ' num2str(i) '=param.FZ(' num2str(i) ',:);' char(10)];
%     new_C_str=[new_F_str 'CX' num2str(i) '=ext_wrenches.CX(' num2str(i) ',:);' 'CY' num2str(i) '=param.CY(' num2str(i) ',:);' 'CZ' num2str(i) '=param.CZ(' num2str(i) ',:);' char(10)];
% 
% end
% 
% 
% %Find and allocate the segments length
% L_ind_first=strfind(fIdstr,'L');
% L_ind_first=L_ind_first(find(L_ind_first<dcl_out_ind ));
% L_ind_first=L_ind_first(find(L_ind_first>dcl_inp_ind ));
% new_L_str=[];
% 
% j=1;
% for i=1:length(L_ind_first)
%     
%     space_next=space_ind(find(space_ind>L_ind_first(i),1));
%     return_next=return_ind(find(return_ind>L_ind_first(i),1));
%         if return_next<space_next
%             new_L_str=[new_L_str fIdstr(L_ind_first(i):return_next-2) '=param.L',fIdstr(L_ind_first(i)+1:return_next-2),';' char(10)];
%         else
%             new_L_str=[new_L_str fIdstr(L_ind_first(i):space_next-1) '=param.L',fIdstr(L_ind_first(i)+1:space_next-1),';' char(10)];
%         end
% end
% 
% 
% 
% 
% 
% %Find and allocate the torque
% %determine which torques should be outputed based on par file()determine
% %where are the active joint
% fileName = ['../Symoro/' modelName,'.par'];
% 
% 
% fId3 = fopen([fileName],'r');
% fId3str = fscanf(fId3,'%c');
% 
% 
% bracket_open_ind=strfind(fId3str,'{');%find all opening and closing brackets containg parameters
% bracket_close_ind=strfind(fId3str,'}');
% Theta_ind=strfind(fId3str,'Theta  = ');
% 
% t_ind=strfind(fId3str(Theta_ind+6:bracket_close_ind(find( bracket_close_ind>Theta_ind,1))),'q')+Theta_ind+5;%find all the t variables (joint angles) contained in the theta vector
% comma_ind=strfind(fId3str(Theta_ind+6:bracket_close_ind(find( bracket_close_ind>Theta_ind,1))),',')+Theta_ind+5;%find all the comma contained in the theta vector
% j=1;
% for i=1:length(comma_ind)-1
%     if isempty(strfind(fId3str(comma_ind(i):comma_ind(i+1)),'q'))==0
%     activejoint(j)=i+1;
%     j=j+1;
%     end
% end
%     
% new_Gamma_str=[];
% for i=1:length(activejoint)
% 
%     new_Gamma_str=[new_Gamma_str 'GAMMA(', num2str(i), ',:)=GAM' num2str(activejoint(i)) ';' char(10)];
% 
% end
% 
% 
% %%Allocate the external wrenche 
% %!!!!! this function is to be applied only for one exit wrenche (at body 0)
% new_EN_str=[];
% 
% new_EN_str=[new_EN_str 'EN(1,:)=E10;EN(2,:)=E20;EN(3,:)=E30;EN(4,:)=N10;EN(5,:)=N20;EN(6,:)=N30;']
% 
% %% Write everything in model_name_dyn2.mat
% fclose(fId);%close original file
% clear fId2
% fId2 = fopen(['../Symoro/' modelName,'_dyn2.m'],'r+')
% 
% fId2str = fscanf(fId2,'%c');
% 
% fId2str = strrep(fId2str,'global','%global');%Comment all global
% 
% 
% 
% fId2str = strrep(fId2str,['% Declaration of %global input variables'],['% Declaration of %global input variables' char(10) new_q_str new_dq_str new_ddq_str new_base0_str new_COM_str ...
%     new_Mass_str new_Inertia_str new_FS_str new_FV_str new_F_str new_C_str]);   
% 
% 
% fId2str = strrep(fId2str,[modelName '_dyn()'],['[GAMMA EN]=' modelName '_dyn2(d,dq,ddq,base0,ext_wrenches,param)']);   
% 
% fId2str = strrep(fId2str,[' ''*'' or ''/'''],[' ''*'' or ''/''' char(10) new_Gamma_str new_EN_str]); 
% 
% 
% fclose(fId2);
% fId2 = fopen(['../Symoro/' modelName,'_dyn2.m'],'w+')
% 
% 
% fprintf(fId2, '%c',[fId2str])
% 
% 
% fclose(fId2);
% 
% 
% fclose(fId3);
% 
% 
% 
% % addpath('Animate_robot_functions\')
% % file=dir('*.dyn');
% %
% %
% % copyfile(file(1).name,[file(1).name 'mat']);%creat a new file specific for matlab
% %
% % %% Replace all comma by semicolon
% % fId = fopen('3D_whole_body_65_frames_dyn.mat','r+')
% % fIdstr = fscanf(fId,'% Declaration ')
% 
% 
% 
% 
