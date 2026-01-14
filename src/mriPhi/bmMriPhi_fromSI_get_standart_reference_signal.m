% Bastien Milani
% CHUV and UNIL
% Lausanne - Switzerland
% May 2023

% argSI must be in the size [nCh, N, nShot] or reshapable.
%
% nSeg must include the SI because the phyical time-list must include the
% SI acquisition. SI's are removed from the lineMask before mitosis.

function [  s_ref, ...
            t_ref, ...
            Fs_ref, ...
            nu_ref, ...
            imNav, ...
            ind_shot_min, ...
            ind_shot_max, ...
            ind_SI_min, ...
            ind_SI_max, ...
            varargout] = bmMriPhi_fromSI_get_standart_reference_signal( rmsSI, ...
                                                                        N, ...
                                                                        nSeg, ...
                                                                        nShot)


% initial -----------------------------------------------------------------
myFigure = figure(  'Name', 'bmMriPhi_refSignal_from_SI', ...
                    'WindowScrollWheelFcn', @myWindowScrollWheelFcn,...
                    'keyreleasefcn', @myKeyReleaseFcn,...
                    'keypressfcn', @myKeyPressFcn,...
                    'WindowButtonDownFcn', @myClickCallback);
                

% initial -----------------------------------------------------------------

rmsSI_mean      = rmsSI(:, :); 
rmsSI_mean      = mean(rmsSI_mean, 1); 
rmsSI_mean      = repmat(rmsSI_mean, [size(rmsSI, 1), 1]); 
rmsSI           = rmsSI./rmsSI_mean; 
% rmsSI           = rmsSI-rmsSI_mean; rmsSI = rmsSI - min(rmsSI(:)); 

control_flag    = false; 
shift_flag      = false; 
s_reverse_flag  = false; 

myColorLimits = [min(rmsSI(:)), max(rmsSI(:))]; 

plot_factor = 1; 

mySize_1    = size(rmsSI, 1);
mySize_2    = size(rmsSI, 2);

s_ref   = [];
t_ref   = [];
Fs_ref  = [];
nu_ref  = [];

ind_shot_min = 1; 
ind_shot_max = nShot; 

ind_SI_min = max(1, fix(N/2) - fix(N/8)); 
ind_SI_max = min(N, fix(N/2) + fix(N/8)); 

ind_imNav_min = 1; 
ind_imNav_max = N; 

s_raw = []; 
update_s_raw(); 
refresh_display;

% END_initial -------------------------------------------------------------

uiwait;


% final -------------------------------------------------------------------

% we work with mirrored-signals
if ~isempty(s_raw)
    
    s_ref = bmMriPhi_fromSI_rawPerShotSignal_to_standartSignal( s_raw, ...
                                                                nSeg, ...
                                                                nShot, ...
                                                                ind_shot_min, ...
                                                                ind_shot_max);     
    
    % We choose time-list with unit step per line
    delta_t     = 1;
    nLine_ref   = size(s_ref, 2);
    t_ref       = (-nLine_ref/2:nLine_ref/2-1)*delta_t;
    T_ref       = nLine_ref*delta_t;
    
    
    delta_nu    = 1/T_ref;
    nu_ref      = (-nLine_ref/2:nLine_ref/2-1)*delta_nu;
    
    
    Fs_ref = bmDFT(s_ref, t_ref, [], 2, 2);

    
    imNav = bmMriPhi_fromSI_imNav(  rmsSI, ...
                                    N,...
                                    nSeg,...
                                    nShot, ...
                                    ind_shot_min,    ind_shot_max, ...
                                    ind_SI_min,      ind_SI_max, ...
                                    ind_imNav_min,   ind_imNav_max, ...
                                    s_reverse_flag   ); 


end



if sum(isnan(s_ref(:))) > 0
    error('NaN detected in the signal');
    return;
end

varargout{1} = ind_SI_min; 
varargout{2} = ind_SI_max; 
varargout{3} = s_reverse_flag; 

return;

% END_final ---------------------------------------------------------------





    function myWindowScrollWheelFcn(src,evnt) % nested function
        if evnt.VerticalScrollCount > 0
            myScrol = max(1, fix(  abs(evnt.VerticalScrollAmount)/  1  ));
            
            delta_ind = ind_SI_max -ind_SI_min; 
            ind_SI_min = ind_SI_min - myScrol;
            ind_SI_max = ind_SI_max - myScrol;
            
            if ind_SI_min < 1
                ind_SI_min = 1;
                ind_SI_max = ind_SI_min + delta_ind; 
            end
        elseif evnt.VerticalScrollCount < 0
            myScrol = max(1, fix(  abs(evnt.VerticalScrollAmount)/  1  ));
            
            delta_ind = ind_SI_max - ind_SI_min; 
            ind_SI_min = ind_SI_min + myScrol;
            ind_SI_max = ind_SI_max + myScrol;
            
            if ind_SI_max > N
                ind_SI_max = N;
                ind_SI_min = ind_SI_max - delta_ind; 
            end
        end
        update_s_raw; 
        refresh_display;
    end



    function myKeyPressFcn(src, command)
        switch lower(command.Key)
            case 'control'     % Ctrl key is pressed
                control_flag = true;
            case 'shift'       % Shift key is pressed
                shift_flag   = true;
            case 'escape'      % Escape key is pressed
                escape_flag  = true;
            case 's'
                s_flag  = true;
            case 'x'
                x_flag  = true;
            case 'n'
                n_flag  = true;
            case 'r'
                if control_flag && ~s_reverse_flag
                    s_reverse_flag  = true;
                    control_flag = false;
                elseif control_flag && s_reverse_flag
                    s_reverse_flag  = false;
                    control_flag = false;
                end
                update_s_raw;
                refresh_display;
                
            case 'downarrow'
                if control_flag
                    plot_factor = plot_factor - 1;
                    if plot_factor < 1
                        plot_factor = 1;
                    end
                else
                    delta_ind = ind_SI_max -ind_SI_min;
                    ind_SI_min = ind_SI_min - 1;
                    ind_SI_max = ind_SI_max - 1;
                    if ind_SI_min < 1
                        ind_SI_min = 1;
                        ind_SI_max = ind_SI_min + delta_ind;
                    end
                end
                update_s_raw;
                refresh_display;
                
            case 'uparrow'
                if control_flag
                    plot_factor = plot_factor + 1; 
                    if plot_factor > N
                       plot_factor = N;  
                    end
                else
                    delta_ind = ind_SI_max -ind_SI_min;
                    ind_SI_min = ind_SI_min + 1;
                    ind_SI_max = ind_SI_max + 1;
                    if ind_SI_max > N
                        ind_SI_max = N;
                        ind_SI_min = ind_SI_max - delta_ind;
                    end
                end
                update_s_raw;
                refresh_display;
                
            case 'leftarrow'
                ind_SI_max = ind_SI_max - 1;
                if ind_SI_max < ind_SI_min
                    ind_SI_max = ind_SI_min;
                end
                update_s_raw;
                refresh_display;
                
            case 'rightarrow'
                ind_SI_max = ind_SI_max + 1;
                if ind_SI_max > N
                    ind_SI_max = N;
                end
                update_s_raw;
                refresh_display;
            case 'e'
                if control_flag
                    imcontrast(myFigure);    % DO NOT REFRESH image after this command
                    control_flag = 0;
                elseif shift_flag
                    myColorLimits = get(gca,'CLim');
                    refresh_display;
                    shift_flag = 0;
                end
        end
    end

    function myKeyReleaseFcn(src, command)
        switch lower(command.Key)
            case 'control'     % Ctrl key is released
                control_flag = false;
            case 'shift'       % Shift key is released
                shift_flag   = false;
            case 'escape'      % Escape key is released
                escape_flag  = false;
            case 's'
                s_flag  = false; 
            case 'x'
                x_flag  = false; 
            case 'n'
                n_flag  = false; 
        end
    end


    function myClickCallback(src, evnt)
        switch get(gcf,'selectiontype')
            case 'alt'% right mouse button click

                p = get_soft_point_from_click;

                [choice_ind, choice_flag] = listdlg('SelectionMode','single', ...
                                                    'ListString', ...
                                                    {'Left Limit', ...
                                                    'Right Limit', ...
                                                    'Bottom Limit', ...
                                                    'Top Limit', ...
                                                    'Bottom View Limit', ...
                                                    'Top View Limit'}); 

                if choice_flag == 0
                    return;
                end

                if choice_ind == 1
                    ind_shot_min = p(2);
                    if ~isempty(ind_shot_max)
                        if ind_shot_min > ind_shot_max
                            ind_shot_min = ind_shot_max;
                        end
                    end

                elseif choice_ind == 2
                    ind_shot_max = p(2);
                    if ~isempty(ind_shot_min)
                        if ind_shot_max < ind_shot_min
                            ind_shot_max = ind_shot_min;
                        end
                    end

                elseif choice_ind == 3
                    ind_SI_min = p(1);
                    if ~isempty(ind_SI_max)
                        if ind_SI_min > ind_SI_max
                            ind_SI_min = ind_SI_max;
                        end
                    end
                    
                elseif choice_ind == 4
                    ind_SI_max = p(1);
                    if ~isempty(ind_SI_min)
                        if ind_SI_max < ind_SI_min
                            ind_SI_max = ind_SI_min;
                        end
                    end

                elseif choice_ind == 5
                    ind_imNav_min = p(1);
                    if ~isempty(ind_imNav_max)
                        if ind_imNav_min > ind_imNav_max
                            ind_imNav_min = ind_imNav_max;
                        end
                    end

                elseif choice_ind == 6
                    ind_imNav_max = p(1);
                    if ~isempty(ind_imNav_min)
                        if ind_imNav_max < ind_imNav_min
                            ind_imNav_max = ind_imNav_min;
                        end
                    end

                end

            case 'normal' % left mouse button click
                1+1; 
            case 'extend'% right mouse button click
                1+1;          
        end
        
        update_s_raw; 
        refresh_display; 
    end

    function update_s_raw()
       
        s = rmsSI(ind_SI_min:ind_SI_max,  :);
        s = mean(s, 1);


        s = s - mean(s(:));
        s = s/std(s(:));
        if s_reverse_flag
            s = -s;
        end
        
        s_raw = s; 

        % [ind_shot_min, ind_shot_max] = bmMriPhi_getIndex_of_firstAndLastExtremum(s_raw); 

    end





    function refresh_display()
        
        figure(myFigure);
        imagesc( rmsSI, [0, 3*mean(rmsSI(:))]  );
        colormap gray
        set(gca, 'YDir', 'normal');
        caxis(myColorLimits)
        hold on

        bmMriPhi_fromSI_plot_signal(s_raw, ind_shot_min, ind_shot_max, ind_SI_min, ind_SI_max, plot_factor); 
        
        if ~isempty(ind_shot_min)
           plot([ind_shot_min, ind_shot_min], [0, mySize_1+1], 'y-') 
        end
        if ~isempty(ind_shot_max)
           plot([ind_shot_max, ind_shot_max], [0, mySize_1+1], 'y-') 
        end
        if ~isempty(ind_SI_min)
           plot([0, mySize_2+1], [ind_SI_min, ind_SI_min], 'y-') 
        end
        if ~isempty(ind_SI_max)
           plot([0, mySize_2+1], [ind_SI_max, ind_SI_max], 'y-')
        end
       
        if ~isempty(ind_imNav_min)
           plot([0, mySize_2+1], [ind_imNav_min, ind_imNav_min], 'g-') 
        end
        if ~isempty(ind_imNav_max)
           plot([0, mySize_2+1], [ind_imNav_max, ind_imNav_max], 'g-')
        end
       
        
        axis([1, mySize_2, 1, mySize_1])
        
        hold off
    end


    function p = get_soft_point_from_click()
        
        myCoordinates = get(gca,'CurrentPoint'); 
        myCoordinates = ceil(myCoordinates(1,1:2)-[0.5 0.5]);     
        
        myX = max(1, myCoordinates(2) );
        myX = min(  mySize_1, myX);
        
        myY = max(1, myCoordinates(1) );
        myY = min(  mySize_2,  myY);
        
        p = [myX, myY]';
        
    end



end


