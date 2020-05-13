function [q, dq, ddq, dddq] = calc_q(splineFit, param)
    switch splineFit.type
        case 'spline'
            q =     fnval(splineFit.sp_qd0, param.t_spline);
            dq =    fnval(splineFit.sp_qd1, param.t_spline);
            ddq =   fnval(splineFit.sp_qd2, param.t_spline);
            dddq =  fnval(splineFit.sp_qd3, param.t_spline);
            
        case 'polyfit'
            dofCount = length(splineFit.sp_qd0);
            tCount = length(param.t_spline);
            q = zeros(dofCount, tCount);
            dq = zeros(dofCount, tCount);
            ddq = zeros(dofCount, tCount);
            dddq = zeros(dofCount, tCount);
            
            for ind_n = 1:length(splineFit.sp_qd0)
                q(ind_n, :) =     polyval(splineFit.sp_qd0{ind_n}, param.t_spline);
                dq(ind_n, :) =    polyval(splineFit.sp_qd1{ind_n}, param.t_spline);
                ddq(ind_n, :) =   polyval(splineFit.sp_qd2{ind_n}, param.t_spline);
                dddq(ind_n, :) =  polyval(splineFit.sp_qd3{ind_n}, param.t_spline);
            end

        case 'piecewise_cds'
            dofCount = length(splineFit.sp_qd0{1}.sp_qd0);
            tCount = length(param.t_spline);
            q = zeros(dofCount, tCount);
            dq = zeros(dofCount, tCount);
            ddq = zeros(dofCount, tCount);
            dddq = zeros(dofCount, tCount);
            
            for ind_n = 1:length(splineFit.sp_qd0) % note that pcds are constructed differently than the other spline types
                param_local = param;
%                 param_local.t_spline = splineFit.sp_qd0{ind_n}.t_spline;
                param_local.t_spline = splineFit.t_spline;
                [q_temp, dq_temp, ddq_temp, dddq_temp] = calc_q(splineFit.sp_qd0{ind_n}, param_local);
                
                switch ind_n
                    case 1
                        ind_start = splineFit.sp_qd0{ind_n}.x_knot(1);
                        ind_end = splineFit.sp_qd0{ind_n}.x_knot(end);
                    otherwise
                        ind_start = splineFit.sp_qd0{ind_n}.x_knot(1)+1; % this keeps it causal
                        ind_end = splineFit.sp_qd0{ind_n}.x_knot(end);
                end
              
                q(:, ind_start:ind_end) = q_temp(:, ind_start:ind_end);
                dq(:, ind_start:ind_end) = dq_temp(:, ind_start:ind_end);
                ddq(:, ind_start:ind_end) = ddq_temp(:, ind_start:ind_end);
                dddq(:, ind_start:ind_end) = dddq_temp(:, ind_start:ind_end);

%                 % plotting purposes
%                 q_temp_holding{ind_n} = q_temp;
%                 dq_temp_holding{ind_n} = dq_temp;
%                 ddq_temp_holding{ind_n} = ddq_temp;
% 
%                 % some symbols for plotting
%  plotSyms = {'x', '^', 'o', '.'};
%  plotColours = {'r', 'b', 'g', 'm', 'k', 'r', 'b', 'g', 'm', 'k', 'r', 'b', 'g', 'm'};
% %                 if 1
% % %                     figure
% %                     subplot(211);
% %                     hold on; plot(q_temp', plotSyms{ind_n});
% %                     plot([ind_end ind_end], ylim, 'k-');
% %                     
% %                     subplot(212);
% %                     hold on; plot(dq_temp', plotSyms{ind_n});
% %                     plot([ind_end ind_end], ylim, 'k-');
% %                 end
            end
            
%             if 0
%                 toplotfeat = q';
%                 toplotfeat2 = q_temp_holding;
%                 
%                 figure;
%                 plot(toplotfeat, 'o'); hold on
%                 yy=axis;
%                 for ind_n=1:length(splineFit.sp_qd0)
%                     ind_start = splineFit.sp_qd0{ind_n}.x_knot(1);
%                     ind_end = splineFit.sp_qd0{ind_n}.x_knot(end);
%                 
%                     plot([ind_start ind_start], yy(3:4));
%                     plot([ind_end ind_end], yy(3:4));
%                     plot(toplotfeat2{ind_n}, plotColours{ind_n});
%                 end
%                 
% %                 for ind_n=1:length(splineFit.sp_qd0)
% %                     plot(splineFit.sp_qd0{ind_n}.x_knot, splineFit.sp_qd0{ind_n}.sp_ed_y, 'x', 'MarkerSize', 20);
% %                 end
%                 
%                 ylim(yy(3:4));
%             end
    end
end