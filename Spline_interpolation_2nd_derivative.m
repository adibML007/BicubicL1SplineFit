function [b,dbdz,caseNum] = Spline_interpolation_2nd_derivative(xinput,zinput)
% THIS IS THE LATEST VERSION.
% This function calculate the b values and db/dz values for
% second-order-derivative based interpolation splines.

% NOTE: the nonuniqueness is solved by choosing the median of the interval
% and a "regTerm". A usual way of solving non-uniqueness is to choose the median of the
% interval and (z3-z1)/(x3-x1). This is the default setting.
% NOTE THAT ALL MIRROR-OPERATION REGARDING REGTERM IS CURRENTLY FOR THE
% DEFAULT SETTING OF REGTERM.

% NOTE: Choice of unique b2 changes values of other b's. Once the choice of
% b2 is changed, check Yu et al. (2010) and Jin et al. (2010) for the formula for other b's.

% NOTE: Choice of subgradient needs further attention and remains not fully
% understood. Old codes sometimes use dz_i and dz_{i+1} interchangeablly
% when dz_i = dz_{i+1}. This issue is important in subgradient search
% method for solving spline fit problems. Considering the fact that subgradient search method
% converges, any choice of subgradient should be OK.

% NOTE: A conjecture: can we approximate any case with "=" by its neighbor
% case? For example, use "+++" to approximate "++=".

% Checked and fixed in 2017
% Last update: 05/04/2017

% Variable definitions
% xinput: horizontal coordinate of data (spline nodes)
% zinput: vertical coordinate of data (spline node values)
% x, z: data in a 5-point window
% b: first-order derivative at spline node
% dbdz: the partial derivatives of b with respect to z
% q0,...,q4: b in a local 5-point window
% g0,...,g4: dbdz in a local 5-point window.

TOL = 1e-7;
CONST00 = ((2 - sqrt(10)) / sqrt(10));  % CONST00 < 0
CONST01 = ((sqrt(10) - 5)/(7 - 2 * sqrt(10)));  % CONST00*CONST01 = 1
CONST02 = ((3 * sqrt(10) - 9)/(7 - 2 * sqrt(10))); % CONST02 > 0; CONST02 = 3/(sqrt(10)+1)
CONST03 = ((2 * sqrt(10) - 2) / (sqrt(10) - 2));
CONST04 = (-3/(7 + sqrt(10)));

n = length(xinput);

b = zeros(1,n);
dbdz = zeros(n,n);

for p = 3:(n-2)
    x = xinput(p-2:p+2);
    z = zinput(p-2:p+2);
    
    h = x(2:end) - x(1:(end-1));
    dz = (z(2:end) - z(1:(end-1)))./h;
    h0 = h(1); h1 = h(2); h2 = h(3); h3 = h(4);
    dz0 = dz(1); dz1 = dz(2); dz2 = dz(3); dz3 = dz(4);
    
    c1 = dz0-dz1;
    c2 = dz3-dz2;
    
    % gradients of "dz". NOTE: mirror operations should involve these
    % subgradients.
    gdz0 = [-1/h0, 1/h0, 0, 0, 0];
    gdz1 = [0, -1/h1, 1/h1, 0, 0];
    gdz2 = [0, 0, -1/h2, 1/h2, 0];
    gdz3 = [0, 0, 0, -1/h3, 1/h3];
    
    % regTerm is the slope between (x1,z1) and (x3,z3). This should not be
    % changed. For other choice of unique b2, define a new variable.
    regTerm = (z(4)-z(2))/(x(4)-x(2));
    gregTerm = [0, -1/(x(4)-x(2)), 0, 1/(x(4)-x(2)), 0];
%     regTerm = 0; % another choice is to choose 0
%     gregTerm = [0 0 0 0 0];
    
    %% case coding
    if dz1 - dz0 < - TOL
        sign1 = 2;
    else if dz1 - dz0 > TOL
            sign1 = 1;
        else
            sign1 = 0;
        end
    end
    
    if dz2 - dz1 < - TOL
        sign2 = 2;
    else if dz2 - dz1 > TOL
            sign2 = 1;
        else
            sign2 = 0;
        end
    end
    
    if dz3 - dz2 < - TOL
        sign3 = 2;
    else if dz3 - dz2 > TOL
            sign3 = 1;
        else
            sign3 = 0;
        end
    end
    
    caseNum = 1 + sign3 + 3 * sign2 + 9 * sign1; % Note the 1+ here.
    
    %% Calculation
    
    switch caseNum
        
        case 1
            % (= = =) unique q2.
            q2 = dz1;
            g2 = gdz1;
            
        case 2
            % (= = +) unique q2. 
            q2 = dz1;
            g2 = gdz1;
            
        case 3
            % (= = -) x-mirror of case 2. x unchanged; z -> -z. dz -> -dz,
            % b -> -b. Replace b_i by -b_i and dz_i by -dz_i. db_idz_j ->
            % db_{i}dz_{j}. gdz_i/d(z_j) -> gdz_{i}/d(z_{j}). So dbdz unchanged; gdz unchanged.
            % More precisely, replace g by -g and gdz by -gdz.
            q2 = dz1;
            g2 = gdz1;
            
        case 4
            % (= + =) nonunique q2. [dz1,dz2].
            q2 = median([regTerm,dz1,dz2]);
            if regTerm < dz1
                g2 = gdz1;
            elseif regTerm > dz2
                g2 = gdz2;
            else
                g2 = gregTerm;
            end
            
        case 5
            % (= + +) unique q2.
            q2 = dz1;
            g2 = gdz1;
            
        case 6
            % (= + -) unique q2.
            q2 = dz1;
            g2 = gdz1;
            
        case 7
            % (= - =) x-mirror of case 4.  dz -> -dz, b -> -b, regTerm ->
            % -regTerm.
            q2 = median([regTerm,dz2,dz1]);
            if regTerm > dz1
                g2 = gdz1;
            elseif regTerm < dz2
                g2 = gdz2;
            else
                g2 = gregTerm;
            end
           
        case 8
            % (= - +) x-mirror of case 6. 
            q2 = dz1;
            g2 = gdz1;
            
        case 9
            % (= - -) x-mirror of case 5
            q2 = dz1;
            g2 = gdz1;
            
        case 10
            % (+ = =) z-mirror of case 2. (x0,x1,x2,x3,x4) ->
            % (-x4,-x3,-x2,-x1,-x0); (z0,z1,z2,z3,z4) ->
            % (z4,z3,z2,z1,z0); (dz0,dz1,dz2,dz3) ->
            % (-dz3,-dz2,-dz1,-dz0); (b0,b1,b2,b3,b4) ->
            % (-b4,-b3,-b2,-b1,-b0). Replace b_{i} by -b_{4-i} and replace
            % dz_{i} by -dz_{3-i}. db_idz_j -> -db_{4-i}dz_{4-j}. gdz_i/d(z_j) -> -gdz_{3-i}/d(z_{4-j}). So replace g_i by
            % g_{4-i} and gdz_i by gdz_{3-i}. The imputation of
            % z_j->z_{4-j} won't impact on g_i because all are linear
            % operations and the imputation will be reversed when mirrored
            % back to g_i.
            q2 = dz2;
            g2 = gdz2;
            
        case 11
            % (+ = +) unique b2.
            q2 = dz1;
            g2 = gdz1;
        
        case 12
            % (+ = -) nonunique b2. b2 in [dz1,dz1+min(CONST00*c1,CONST00*c2)].
            lb = dz1; ub = dz1+min(CONST00*c1,CONST00*c2);
            q2 = median([regTerm,lb,ub]);
            if regTerm <= lb
                g2 = gdz1;
            elseif regTerm < ub
                g2 = gregTerm;
            else
                if CONST00*c1 <= CONST00*c2
                    g2 = gdz1 + CONST00 * (gdz0 - gdz1);
                else
                    g2 = gdz1 + CONST00 * (gdz3 - gdz2);
                end
            end
            
        case 13
            % (+ + =) z-mirror of case 5. Replace b_{i} by -b_{4-i} and replace
            % dz_{i} by -dz_{3-i}.
            q2 = dz2;
            g2 = gdz2;
           
        case 14
            % (+ + +) 4 subcases
            if (dz3 - dz0) >= CONST03 * (dz2 - dz1)  % subcase 14-1. This condition is equivalent to dz2 - dz1 <= -CONST00*(-c1+c2)
                % nonunique b2.
                % [max(dz1,dz2+CONST00*c2),min(dz1+CONST00*c1,dz2)]. 
                lb = dz2+CONST00*(dz3-dz2); ub = dz1+CONST00*(dz0-dz1);
                glb = gdz2+CONST00*(gdz3-gdz2); gub = gdz1+CONST00*(gdz0-gdz1);
                if regTerm < max(dz1,lb)
                    q2 = max(dz1,lb);
                    if dz1 <= lb
                        g2 = glb;
                    else
                        g2 = gdz1;
                    end
                elseif regTerm > min(dz2,ub)
                    q2 = min(dz2,ub);
                    if dz2 <= ub
                        g2 = gdz2;
                    else
                        g2 = gub;
                    end
                else
                    q2 = regTerm;
                    g2 = gregTerm;
                end
            else % subcases 14-2, 14-3, 14-4
                q2Bound1 = dz1 - (dz0 - dz1) / 2;
                q2Bound2 = dz1 - (dz0 - dz1) * 2;
                gq2Bound1 = gdz1 - (gdz0 - gdz1) / 2;
                gq2Bound2 = gdz1 - (gdz0 - gdz1) * 2;
                
                if -2 * (q2Bound1 - dz2) < (dz3 - dz2) % subcase 14-2. unique b2; analytic (line search)
                    alpha = (dz3 - dz2) / (dz3 - dz2 + dz1 - dz0);
                    galpha = ((gdz3 - gdz2)*(dz3 - dz2 + dz1 - dz0)- (gdz3 - gdz2 + gdz1 - gdz0)*(dz3 - dz2))/(dz3 - dz2 + dz1 - dz0)^2 ;
                    q2 = alpha * dz1 + (1.0 - alpha) * dz2;
                    g2 = galpha * dz1 + alpha * gdz1 - galpha * dz2 + (1-alpha) * gdz2;
                elseif dz3 - dz2 < -(q2Bound2 - dz2) / 2.0 % subcase 14-4. unique b2; analytic (line search)
                    alpha = (dz3 - dz2) / (dz3 - dz2 + dz1 - dz0);
                    galpha = ((gdz3 - gdz2)*(dz3 - dz2 + dz1 - dz0)- (gdz3 - gdz2 + gdz1 - gdz0)*(dz3 - dz2))/(dz3 - dz2 + dz1 - dz0)^2 ;
                    q2 = alpha * dz1 + (1.0 - alpha) * dz2;
                    g2 = galpha * dz1 + alpha * gdz1 - galpha * dz2 + (1-alpha) * gdz2;
                else % subcase 14-3. nonunique b2. [max(dz1-c1/2,dz2-2*c2),min(dz1-2*c1,dz2-c2/2)]
                    lb = dz2 - 2 * (dz3 - dz2);
                    ub = dz2 - (dz3 - dz2) / 2;
                    glb = gdz2 - 2*(gdz3-gdz2);
                    gub = gdz2 - (gdz3 - gdz2)/2;
                    if regTerm < max(q2Bound1, lb)
                        q2 = max(q2Bound1, lb);
                        if q2Bound1 < lb
                            g2 = glb;
                        else
                            g2 = gq2Bound1;
                        end
                    elseif regTerm > min(q2Bound2, ub)
                        q2 = min(q2Bound2, ub);
                        if q2Bound2 < ub
                            g2 = gq2Bound2;
                        else
                            g2 = gub;
                        end
                    else
                        q2 = regTerm;
                        g2 = gregTerm;
                    end
                end
            end
            
        case 15
            % (+ + -) 3 subcases
            if dz0 - dz1 < CONST04 * (dz2 - dz1) % subcases 15-1 and 15-2
                if dz2 - dz1 <= CONST00 * (dz0-dz1)  % subcase 15-1, nonunique b2. [dz2,min(dz1+CONST00*c1,dz2+CONST00*c2)].
                    ub = min(dz1+CONST00*(dz0-dz1),dz2+CONST00*(dz3-dz2));
                    q2 = median([regTerm,dz2,ub]);
                    if regTerm <= dz2
                        g2 = gdz2;
                    elseif regTerm <= ub
                        g2 = gregTerm;
                    else
                        if dz1+CONST00*(dz0-dz1) <= dz2+CONST00*(dz3-dz2)
                            g2 = gdz1 + CONST00 * (gdz0 - gdz1);
                        else
                            g2 = gdz2 + CONST00 * (gdz3 - gdz2);
                        end
                    end
                else % subcase 15-2, unique b2
                    q2 = dz2;
                    g2 = gdz2;
                end
            else  % subcase 15-3; (dz0 - dz1)/CONST04 < dz2-dz1
                q2Bound1 = dz1 + (dz0 - dz1) / CONST04;
                gq2Bound1 = gdz1 + (gdz0 - gdz1)/CONST04;
                if dz3 - dz2 < CONST02 * (q2Bound1 - dz2) % subcase 15-3-1, unique b2. See Qingwei's handwritten doc.
                    q2 = q2Bound1;
                    g2 = gq2Bound1;
                else % subcase 15-3-2, analytic (line search). See Qingwei's handwritten doc.
                    alpha = (dz3 - dz2) / (dz3 - dz2 - dz1 + dz0);
                    galpha = ( (gdz3 - gdz2)*(dz3 - dz2 - dz1 + dz0) - (gdz3-gdz2-gdz1+gdz0)*(dz3 - dz2) )/(dz3 - dz2 - dz1 + dz0)^2;
                    q2 = alpha * (2 * dz1 - dz0) + (1 - alpha) * (2 * dz2 - dz3);
                    g2 = galpha * (2*dz1-dz0) + alpha*(2*gdz1-gdz0) - galpha*(2.0 * dz2 - dz3) + (1.0 - alpha) * (2.0 * gdz2 - gdz3);
                end
            end
            
        case 16
            % (+ - =). (0,0)-rotation of case 6 (equivalent to x-mirror
            % then z-mirror). (x0,x1,x2,x3,x4) -> (-x4,-x3,-x2,-x1,-x0);
            % (z0,z1,z2,z3,z4) -> (-z4,-z3,-z2,-z1,-z0); (dz0,dz1,dz2,dz3)
            % -> (dz3,dz2,dz1,dz0); regTerm -> regTerm. (b0,b1,b2,b3,b4) ->
            % (b4,b3,b2,b1,b0). Replace b_{i} by b_{4-i}. Replace dz_i by
            % dz_{3-i}. db_idz_j -> -db_{4-i}dz_{4-j}. gdz_i/d(z_j) -> -gdz_{3-i}/d(z_{4-j}).
            q2 = dz2;
            g2 = gdz2;
            
        case 17
            % (+ - +) 2 subcases.
            if dz1 - dz0 + dz3 - dz2 >= CONST02 * (dz1 - dz2) % subcase 17-2 nonunique b2.
                % optimal interval of b2 is [max(dz2, dz1 + (dz0 - dz1) / CONST02),min(dz1, dz2 + (dz3 - dz2) / CONST02)]
                q2Bound1 = dz1 + (dz0 - dz1) / CONST02;
                q2Bound2 = dz2 + (dz3 - dz2) / CONST02;
                gq2Bound1 = gdz1 + (gdz0 - gdz1) / CONST02;
                gq2Bound2 = gdz2 + (gdz3 - gdz2) / CONST02;
                if regTerm < max(dz2,q2Bound1)
                    q2 = max(dz2,q2Bound1);
                    if dz2 <= q2Bound1
                        g2 = gq2Bound1;
                    else
                        g2 = gdz2;
                    end
                elseif regTerm > min(dz1,q2Bound2)
                    q2 = min(dz1,q2Bound2);
                    if dz1 <= q2Bound2
                        g2 = gdz1;
                    else
                        g2 = gq2Bound2;
                    end
                else
                    q2 = regTerm;
                    g2 = gregTerm;
                end
            else % subcase 17-1, unique b2
                alpha = (dz3 - dz2) / (dz3 - dz2 + dz1 - dz0);
                galpha = ( (gdz3 - gdz2)*(dz3 - dz2 + dz1 - dz0)-(gdz3 - gdz2 + gdz1 - gdz0)*(dz3 - dz2) ) / (dz3 - dz2 + dz1 - dz0)^2;
                q2 = alpha * dz1 + (1.0 - alpha) * dz2;
                g2 = galpha * dz1 + alpha*gdz1 - galpha*dz2 + (1-alpha)*gdz2;
            end
            
        case 18
            % (+ - -) (0,0)-rotation of case 15 (equivalent to x-mirror
            % then z-mirror). (x0,x1,x2,x3,x4) -> (-x4,-x3,-x2,-x1,-x0);
            % (z0,z1,z2,z3,z4) -> (-z4,-z3,-z2,-z1,-z0); (dz0,dz1,dz2,dz3)
            % -> (dz3,dz2,dz1,dz0); regTerm -> regTerm. (b0,b1,b2,b3,b4) ->
            % (b4,b3,b2,b1,b0). Replace q_{i} by q_{4-i}. Replace dz_i by
            % dz_{3-i}. Replace g_i by -g_{4-i} and gdz_i by -gdz_{3-i}.
            dz0Temp = -dz0; gdz0Temp = gdz0;
            dz1Temp = -dz1; gdz1Temp = gdz1;
            dz2Temp = -dz2; gdz2Temp = gdz2;
            dz3Temp = -dz3; gdz3Temp = gdz3;
            dz0 = -dz3Temp; gdz0 = -gdz3Temp;
            dz1 = -dz2Temp; gdz1 = -gdz2Temp;
            dz2 = -dz1Temp; gdz2 = -gdz1Temp;
            dz3 = -dz0Temp; gdz3 = -gdz0Temp;
            % regTerm no change
            gregTerm = -gregTerm; % gregTerm is changed.
            % Copy of case 15
            if dz0 - dz1 < CONST04 * (dz2 - dz1) % subcases 15-1 and 15-2
                if dz2 - dz1 <= CONST00 * (dz0-dz1)  % subcase 15-1, nonunique b2. [dz2,min(dz1+CONST00*c1,dz2+CONST00*c2)].
                    ub = min(dz1+CONST00*(dz0-dz1),dz2+CONST00*(dz3-dz2));
                    q2 = median([regTerm,dz2,ub]);
                    if regTerm <= dz2
                        g2 = gdz2;
                    elseif regTerm <= ub
                        g2 = gregTerm;
                    else
                        if dz1+CONST00*(dz0-dz1) <= dz2+CONST00*(dz3-dz2)
                            g2 = gdz1 + CONST00 * (gdz0 - gdz1);
                        else
                            g2 = gdz2 + CONST00 * (gdz3 - gdz2);
                        end
                    end
                else % subcase 15-2, unique b2
                    q2 = dz2;
                    g2 = gdz2;
                end
            else  % subcase 15-3; (dz0 - dz1)/CONST04 < dz2-dz1
                q2Bound1 = dz1 + (dz0 - dz1) / CONST04;
                gq2Bound1 = gdz1 + (gdz0 - gdz1)/CONST04;
                if dz3 - dz2 < CONST02 * (q2Bound1 - dz2) % subcase 15-3-1, unique b2. See Qingwei's handwritten doc.
                    q2 = q2Bound1;
                    g2 = gq2Bound1;
                else % subcase 15-3-2, analytic (line search). See Qingwei's handwritten doc.
                    alpha = (dz3 - dz2) / (dz3 - dz2 - dz1 + dz0);
                    galpha = ( (gdz3 - gdz2)*(dz3 - dz2 - dz1 + dz0) - (gdz3-gdz2-gdz1+gdz0)*(dz3 - dz2) )/(dz3 - dz2 - dz1 + dz0)^2;
                    q2 = alpha * (2 * dz1 - dz0) + (1 - alpha) * (2 * dz2 - dz3);
                    g2 = galpha * (2*dz1-dz0) + alpha*(2*gdz1-gdz0) - galpha*(2.0 * dz2 - dz3) + (1.0 - alpha) * (2.0 * gdz2 - gdz3);
                end
            end
            % Rotation back. q2 = q2 and g2 = -g2.
            g2 = -g2;
            
        case 19
            % (- = =) (0,0)-rotation of case 2. Replace q_{i} by q_{4-i}. Replace dz_i by
            % dz_{3-i}. unique b2.
            q2 = dz2;
            g2 = gdz2;
            
        case 20
            % (- = +) x-mirror of case 12. dz -> -dz, b -> -b, regTerm ->
            % -regTerm. Nonunique b2. q2 =
            % [dz1+max(CONST00*(c2),CONST00*(c1)),dz1]. It can also be
            % viewed as z-mirror of case 12 and the final value of q2 will
            % be the same but g2 will not be. (subgradient issue)
            dz0 = -dz0; dz1 = -dz1; dz2 = -dz2; dz3 = -dz3;
            c1 = -c1; c2 = -c2; regTerm = -regTerm;
            % copy case 12
            lb = dz1; ub = dz1+min(CONST00*c1,CONST00*c2);
            q2 = median([regTerm,lb,ub]);
            if regTerm <= lb
                g2 = gdz1;
            elseif regTerm < ub
                g2 = gregTerm;
            else
                if CONST00*c1 <= CONST00*c2
                    g2 = gdz1 + CONST00 * (gdz0 - gdz1);
                else
                    g2 = gdz1 + CONST00 * (gdz3 - gdz2);
                end
            end
            % mirror back, q2 -> -q2 but g2 is unchanged.
            q2 = -q2;
            
        case 21
            % (- = -) x-mirror of case 11. unique b2.
            q2 = dz1;
            g2 = gdz1;
            
        case 22
            % (- + =) z-mirror of case 6. unique b2.
            q2 = dz2;
            g2 = gdz2;
           
        case 23
            % (- + +) z-mirror of case 15
            dz0Temp = dz0; gdz0Temp = gdz0;
            dz1Temp = dz1; gdz1Temp = gdz1;
            dz2Temp = dz2; gdz2Temp = gdz2;
            dz3Temp = dz3; gdz3Temp = gdz3;
            dz0 = -dz3Temp; gdz0 = -gdz3Temp;
            dz1 = -dz2Temp; gdz1 = -gdz2Temp;
            dz2 = -dz1Temp; gdz2 = -gdz1Temp;
            dz3 = -dz0Temp; gdz3 = -gdz0Temp;
            regTerm = - regTerm; gregTerm = -gregTerm;
            % copy of case 15
            if dz0 - dz1 < CONST04 * (dz2 - dz1) % subcases 15-1 and 15-2
                if dz2 - dz1 <= CONST00 * (dz0-dz1)  % subcase 15-1, nonunique b2. [dz2,min(dz1+CONST00*c1,dz2+CONST00*c2)].
                    ub = min(dz1+CONST00*(dz0-dz1),dz2+CONST00*(dz3-dz2));
                    q2 = median([regTerm,dz2,ub]);
                    if regTerm <= dz2
                        g2 = gdz2;
                    elseif regTerm <= ub
                        g2 = gregTerm;
                    else
                        if dz1+CONST00*(dz0-dz1) <= dz2+CONST00*(dz3-dz2)
                            g2 = gdz1 + CONST00 * (gdz0 - gdz1);
                        else
                            g2 = gdz2 + CONST00 * (gdz3 - gdz2);
                        end
                    end
                else % subcase 15-2, unique b2
                    q2 = dz2;
                    g2 = gdz2;
                end
            else  % subcase 15-3; (dz0 - dz1)/CONST04 < dz2-dz1
                q2Bound1 = dz1 + (dz0 - dz1) / CONST04;
                gq2Bound1 = gdz1 + (gdz0 - gdz1)/CONST04;
                if dz3 - dz2 < CONST02 * (q2Bound1 - dz2) % subcase 15-3-1, unique b2. See Qingwei's handwritten doc.
                    q2 = q2Bound1;
                    g2 = gq2Bound1;
                else % subcase 15-3-2, analytic (line search). See Qingwei's handwritten doc.
                    alpha = (dz3 - dz2) / (dz3 - dz2 - dz1 + dz0);
                    galpha = ( (gdz3 - gdz2)*(dz3 - dz2 - dz1 + dz0) - (gdz3-gdz2-gdz1+gdz0)*(dz3 - dz2) )/(dz3 - dz2 - dz1 + dz0)^2;
                    q2 = alpha * (2 * dz1 - dz0) + (1 - alpha) * (2 * dz2 - dz3);
                    g2 = galpha * (2*dz1-dz0) + alpha*(2*gdz1-gdz0) - galpha*(2.0 * dz2 - dz3) + (1.0 - alpha) * (2.0 * gdz2 - gdz3);
                end
            end
            % mirror back
            q2 = -q2;
            g2 = -g2;
            
        case 24
            % (- + -). x-mirror of case 17
            dz0 = -dz0;
            dz1 = -dz1;
            dz2 = -dz2;
            dz3 = -dz3;
            regTerm = -regTerm;
            % copy of case 17
            if dz1 - dz0 + dz3 - dz2 >= CONST02 * (dz1 - dz2) % subcase 17-2 nonunique b2.
                % optimal interval of b2 is [max(dz2, dz1 + (dz0 - dz1) / CONST02),min(dz1, dz2 + (dz3 - dz2) / CONST02)]
                q2Bound1 = dz1 + (dz0 - dz1) / CONST02;
                q2Bound2 = dz2 + (dz3 - dz2) / CONST02;
                gq2Bound1 = gdz1 + (gdz0 - gdz1) / CONST02;
                gq2Bound2 = gdz2 + (gdz3 - gdz2) / CONST02;
                if regTerm < max(dz2,q2Bound1)
                    q2 = max(dz2,q2Bound1);
                    if dz2 <= q2Bound1
                        g2 = gq2Bound1;
                    else
                        g2 = gdz2;
                    end
                elseif regTerm > min(dz1,q2Bound2)
                    q2 = min(dz1,q2Bound2);
                    if dz1 <= q2Bound2
                        g2 = gdz1;
                    else
                        g2 = gq2Bound2;
                    end
                else
                    q2 = regTerm;
                    g2 = gregTerm;
                end
            else % subcase 17-1, unique b2
                alpha = (dz3 - dz2) / (dz3 - dz2 + dz1 - dz0);
                galpha = ( (gdz3 - gdz2)*(dz3 - dz2 + dz1 - dz0)-(gdz3 - gdz2 + gdz1 - gdz0)*(dz3 - dz2) ) / (dz3 - dz2 + dz1 - dz0)^2;
                q2 = alpha * dz1 + (1.0 - alpha) * dz2;
                g2 = galpha * dz1 + alpha*gdz1 - galpha*dz2 + (1-alpha)*gdz2;
            end
            % mirror back
            q2 = -q2;
            
        case 25
            % (- - =). (0,0)-rotation of case 5. Replace q_{i} by q_{4-i}. Replace dz_i by
            % dz_{3-i}. unique b2.
            q2 = dz2;
            g2 = gdz2;
                        
        case 26
            % (- - +). x-mirror of case 15.
            dz0 = -dz0;
            dz1 = -dz1;
            dz2 = -dz2;
            dz3 = -dz3;
            regTerm = -regTerm;
            % copy of case 15
            if dz0 - dz1 < CONST04 * (dz2 - dz1) % subcases 15-1 and 15-2
                if dz2 - dz1 <= CONST00 * (dz0-dz1)  % subcase 15-1, nonunique b2. [dz2,min(dz1+CONST00*c1,dz2+CONST00*c2)].
                    ub = min(dz1+CONST00*(dz0-dz1),dz2+CONST00*(dz3-dz2));
                    q2 = median([regTerm,dz2,ub]);
                    if regTerm <= dz2
                        g2 = gdz2;
                    elseif regTerm <= ub
                        g2 = gregTerm;
                    else
                        if dz1+CONST00*(dz0-dz1) <= dz2+CONST00*(dz3-dz2)
                            g2 = gdz1 + CONST00 * (gdz0 - gdz1);
                        else
                            g2 = gdz2 + CONST00 * (gdz3 - gdz2);
                        end
                    end
                else % subcase 15-2, unique b2
                    q2 = dz2;
                    g2 = gdz2;
                end
            else  % subcase 15-3; (dz0 - dz1)/CONST04 < dz2-dz1
                q2Bound1 = dz1 + (dz0 - dz1) / CONST04;
                gq2Bound1 = gdz1 + (gdz0 - gdz1)/CONST04;
                if dz3 - dz2 < CONST02 * (q2Bound1 - dz2) % subcase 15-3-1, unique b2. See Qingwei's handwritten doc.
                    q2 = q2Bound1;
                    g2 = gq2Bound1;
                else % subcase 15-3-2, analytic (line search). See Qingwei's handwritten doc.
                    alpha = (dz3 - dz2) / (dz3 - dz2 - dz1 + dz0);
                    galpha = ( (gdz3 - gdz2)*(dz3 - dz2 - dz1 + dz0) - (gdz3-gdz2-gdz1+gdz0)*(dz3 - dz2) )/(dz3 - dz2 - dz1 + dz0)^2;
                    q2 = alpha * (2 * dz1 - dz0) + (1 - alpha) * (2 * dz2 - dz3);
                    g2 = galpha * (2*dz1-dz0) + alpha*(2*gdz1-gdz0) - galpha*(2.0 * dz2 - dz3) + (1.0 - alpha) * (2.0 * gdz2 - gdz3);
                end
            end
            % mirror back
            q2 = -q2;
                        
        case 27
            % (- - -). x-mirror of case 14
            dz0 = - dz0;
            dz1 = - dz1;
            dz2 = - dz2;
            dz3 = - dz3;
            regTerm = - regTerm;
            % copy of case 14
            if (dz3 - dz0) >= CONST03 * (dz2 - dz1)  % subcase 14-1. This condition is equivalent to dz2 - dz1 <= -CONST00*(-c1+c2)
                % nonunique b2.
                % [max(dz1,dz2+CONST00*c2),min(dz1+CONST00*c1,dz2)]. 
                lb = dz2+CONST00*(dz3-dz2); ub = dz1+CONST00*(dz0-dz1);
                glb = gdz2+CONST00*(gdz3-gdz2); gub = gdz1+CONST00*(gdz0-gdz1);
                if regTerm < max(dz1,lb)
                    q2 = max(dz1,lb);
                    if dz1 <= lb
                        g2 = glb;
                    else
                        g2 = gdz1;
                    end
                elseif regTerm > min(dz2,ub)
                    q2 = min(dz2,ub);
                    if dz2 <= ub
                        g2 = gdz2;
                    else
                        g2 = gub;
                    end
                else
                    q2 = regTerm;
                    g2 = gregTerm;
                end
            else % subcases 14-2, 14-3, 14-4
                q2Bound1 = dz1 - (dz0 - dz1) / 2;
                q2Bound2 = dz1 - (dz0 - dz1) * 2;
                gq2Bound1 = gdz1 - (gdz0 - gdz1) / 2;
                gq2Bound2 = gdz1 - (gdz0 - gdz1) * 2;
                
                if -2 * (q2Bound1 - dz2) < (dz3 - dz2) % subcase 14-2. unique b2; analytic (line search)
                    alpha = (dz3 - dz2) / (dz3 - dz2 + dz1 - dz0);
                    galpha = ((gdz3 - gdz2)*(dz3 - dz2 + dz1 - dz0)- (gdz3 - gdz2 + gdz1 - gdz0)*(dz3 - dz2))/(dz3 - dz2 + dz1 - dz0)^2 ;
                    q2 = alpha * dz1 + (1.0 - alpha) * dz2;
                    g2 = galpha * dz1 + alpha * gdz1 - galpha * dz2 + (1-alpha) * gdz2;
                elseif dz3 - dz2 < -(q2Bound2 - dz2) / 2.0 % subcase 14-4. unique b2; analytic (line search)
                    alpha = (dz3 - dz2) / (dz3 - dz2 + dz1 - dz0);
                    galpha = ((gdz3 - gdz2)*(dz3 - dz2 + dz1 - dz0)- (gdz3 - gdz2 + gdz1 - gdz0)*(dz3 - dz2))/(dz3 - dz2 + dz1 - dz0)^2 ;
                    q2 = alpha * dz1 + (1.0 - alpha) * dz2;
                    g2 = galpha * dz1 + alpha * gdz1 - galpha * dz2 + (1-alpha) * gdz2;
                else % subcase 14-3. nonunique b2. [max(dz1-c1/2,dz2-2*c2),min(dz1-2*c1,dz2-c2/2)]
                    lb = dz2 - 2 * (dz3 - dz2);
                    ub = dz2 - (dz3 - dz2) / 2;
                    glb = gdz2 - 2*(gdz3-gdz2);
                    gub = gdz2 - (gdz3 - gdz2)/2;
                    if regTerm < max(q2Bound1, lb)
                        q2 = max(q2Bound1, lb);
                        if q2Bound1 < lb
                            g2 = glb;
                        else
                            g2 = gq2Bound1;
                        end
                    elseif regTerm > min(q2Bound2, ub)
                        q2 = min(q2Bound2, ub);
                        if q2Bound2 < ub
                            g2 = gq2Bound2;
                        else
                            g2 = gub;
                        end
                    else
                        q2 = regTerm;
                        g2 = gregTerm;
                    end
                end
            end
            % mirror back
            q2 = -q2;
            
    end
    % calculate b1
    % restore dz_i
    dz0 = dz(1); dz1 = dz(2); dz2 = dz(3); dz3 = dz(4);
    q1 = dz1 + median([CONST01*(q2-dz1),CONST02*(q2-dz1),dz0-dz1]);
    if CONST01*(q2-dz1) <= CONST02*(q2-dz1)
        if CONST02*(q2-dz1) <= dz0-dz1
            g1 = gdz1 + CONST02*(g2-gdz1);
        else
            if CONST01*(q2-dz1) <= dz0-dz1
                g1 = gdz0;
            else
                g1 = gdz1 + CONST01*(g2-gdz1);
            end
        end
    else
        if CONST01*(q2-dz1) <= dz0-dz1
            g1 = gdz1 + CONST01*(g2-gdz1);
        else
            if CONST02*(q2-dz1) <= dz0-dz1
                g1 = gdz0;
            else
                g1 = gdz1 + CONST02*(g2-gdz1);
            end
        end
    end
    % calculate b3
    q3 = dz2 + median([CONST01*(q2-dz2),CONST02*(q2-dz2),dz3-dz2]);
    if CONST01*(q2-dz2) <= CONST02*(q2-dz2)
        if CONST02*(q2-dz2) <= dz3-dz2
            g3 = gdz2 + CONST02*(g2-gdz2);
        else
            if CONST01*(q2-dz2) <= dz3-dz2
                g3 = gdz3;
            else
                g3 = gdz2 + CONST01*(g2-gdz2);
            end
        end
    else
        if CONST01*(q2-dz2) <= dz3-dz2
            g3 = gdz2 + CONST01*(g2-gdz2);
        else
            if CONST02*(q2-dz2) <= dz3-dz2
                g3 = gdz3;
            else
                g3 = gdz2 + CONST02*(g2-gdz2);
            end
        end
    end
    % calculate b0 and b4
    q0 = dz0 + CONST00*(q1-dz0);
    q4 = dz3 + CONST00*(q3-dz3);
    g0 = gdz0 + CONST00*(g1-gdz0);
    g4 = gdz3 + CONST00*(g3-gdz3);
    
    if p == 3
        b(1) = q0;
        b(2) = q1;
        dbdz(1,p-2:p+2) = g0;
        dbdz(2,p-2:p+2) = g1;
    end
    if p == n-2
        b(n-1) = q3;
        b(n) = q4;
        dbdz(n-1,p-2:p+2) = g3;
        dbdz(n,p-2:p+2) = g4;
    end
    b(p) = q2;
    dbdz(p,p-2:p+2) = g2;
end