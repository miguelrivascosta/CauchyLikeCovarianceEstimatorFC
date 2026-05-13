classdef SheinvaldMethod
    properties
        n
		nrf
        d
        lambda
        nrf2
		Gv
		Q
		M
        L
        Av_grid
        theta_grid_deg;
        epsilon
    end
    methods
        function o = SheinvaldMethod(d,lambda,n,nrf,L,theta0,theta1,theta_step)
			o.n = n;
			o.nrf = nrf;
            o.nrf2 = nrf*nrf;
            o.L = L;
            o.d = d;
            o.lambda = lambda;
            [o.Gv,o.M] = o.gen_codebook();
            o.theta_grid_deg = theta0:theta_step:theta1;
            o.Q = length(o.theta_grid_deg);
            o.epsilon = 0.001;
            Av_grid = zeros(o.nrf,o.Q,o.M);
            for i = 1 : o.M
                Av_grid(:,:,i) = o.Gv(:,:,i)*o.a_ula(2*pi*d/lambda*sind(o.theta_grid_deg));
            end
            o.Av_grid = Av_grid;
        end
		


        function [Gv,M] = gen_codebook(o)
            psi_u = 2*pi/o.n*((0:o.n-1)-o.n/2);
            B = 1/sqrt(o.n)*ToolsObj.a_ula(o.n,psi_u);
            
            if o.nrf==o.n
	            M = 1;
            else
	            M = ceil(o.n/(o.nrf-1));
            end
            
            IDX = zeros(M,o.nrf);
            IDX(1,:) = 0:o.nrf-1;
            for i = 1 : M-1
	            IDX(i+1,:) = mod(IDX(i,:) + o.nrf-1,o.n);
            end
            
            Gv = zeros(o.nrf,o.n,M);
            for i = 1 : M
	            Gv(:,:,i) = B(:,IDX(i,:)+1)'; 
            end
        end

		

		function y = a_ula(o,psi)
			y = exp(1j*(0:o.n-1).'.*psi);
        end

        function theta_hat_deg = estimate_theta_deg(o,X,Km)
            vec_I = repmat(sqrt(Km)*reshape(eye(o.nrf),[o.nrf2,1]),[o.M,1]);
            I = eye(o.M*o.nrf2);
            Xm = reshape(X,[o.n,Km,o.M]);
		    Ym = pagemtimes(o.Gv,Xm);
		    Rv_hat = pagemtimes(Ym,pagectranspose(Ym))/Km;
		    
            Tv = zeros(o.nrf,o.nrf,o.M);
            inv_Rv_hat = pageinv(Rv_hat);
            for v = 1 : o.M
                Tv(:,:,v) = nthroot(Km,4)*sqrtm(inv_Rv_hat(:,:,v));
            end
    
		    rhov_vec = sqrt(Km)*reshape(inv_Rv_hat,[o.M*o.nrf2,1]); 
            Av_tilde = pagemtimes(Tv,o.Av_grid);
		    
            loss = zeros(1,o.Q);
            
            
            idx_signals = zeros(1,o.L);
            A_cal_mat = zeros(o.M*o.nrf2,o.Q);
    
            for q = 1 : o.Q
                A_cal = zeros(o.nrf2,o.M);
                for v = 1 : o.M
                    A_cal(:,v) = kron(conj(Av_tilde(:,q,v)),Av_tilde(:,q,v));
                end
                A_cal_mat(:,q) = A_cal(:);
    
                A_calq = [A_cal_mat(:,q),rhov_vec];
                A_cal_tildeq = [real(A_calq );imag(A_calq )];
			    Ptheta = I-A_calq*((A_cal_tildeq.'*A_cal_tildeq)\real(A_calq).');
                loss(q) = norm(Ptheta*vec_I);
            end
            [~,idx_signals(1)] = min(loss);
            
            if o.L>1
                for l=2:o.L
                    loss = zeros(1,o.Q);
                    fixed_signals = idx_signals(1:l-1);
                    A_cal = zeros(o.M*o.nrf2,l+1);
                    A_cal(:,1:l-1) = A_cal_mat(:,fixed_signals);
                    A_cal(:,end) = rhov_vec;
                    for q = 1 : o.Q
                        if ismember(q,fixed_signals)
                            loss(q)=NaN;
                            continue
                        end
                        A_cal(:,l) = A_cal_mat(:,q);
                        A_cal_tilde = [real(A_cal);imag(A_cal)];
			            Ptheta = I-A_cal*((A_cal_tilde.'*A_cal_tilde)\real(A_cal).');
                        loss(q) = norm(Ptheta*vec_I);
                    end
                    [~,idx_signals(l)] = min(loss);
                end
                
                theta_hat_old = o.theta_grid_deg(sort(idx_signals,'ascend'));
    
                l = 1;
                A_cal = [A_cal_mat(:,idx_signals),rhov_vec];
                while(1)
                    fixed_signals = idx_signals(1:o.L~=l);
                    loss = zeros(1,o.Q);
                    for q = 1 : o.Q
                        if ismember(q,fixed_signals)
                            loss(q)=NaN;
                            continue
                        end
                        A_cal(:,l) = A_cal_mat(:,q);
                        A_cal_tilde = [real(A_cal);imag(A_cal)];
		                Ptheta = I-A_cal*((A_cal_tilde.'*A_cal_tilde)\real(A_cal).');
                        loss(q) = norm(Ptheta*vec_I);
                    end
                    [~,new_idx] = min(loss);
                    idx_signals(l) = new_idx;


                    A_cal = [A_cal_mat(:,idx_signals),rhov_vec];
                    theta_hat = o.theta_grid_deg(sort(idx_signals,'ascend'));

                    if (abs(theta_hat_old-theta_hat)<=o.epsilon)
                        break
                    end

                    theta_hat_old = theta_hat;
                    l = l + 1;
                    if l>o.L
                        l = 1;
                    end
                end
            end
            theta_hat_deg = o.theta_grid_deg(sort(idx_signals,'ascend'));
        end

    end
end
