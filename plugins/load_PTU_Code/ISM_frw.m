function W_ISM_img = ISM_frw(data, psf, eps)    
    
    if nargin < 3
        %regulation parameter to prevent divergence
        eps = 0.07;
    end
    %iamge size 
    [nx,ny] = size(data,[1 2]);
    if nx<1000 || ny<1000
        Nx = 1000;
        Ny = 1000;
    else
        Nx = nx;
        Ny = ny;
    end

    %offset to zero
    psf = psf(1:2:end,1:2:end);
    psf = psf - min(psf,[],'all');
    %remove possible NaNs
    psf(isnan(psf)) = 0.0;
    %Fouier transform
    Fpsf = abs(fftshift(fft2(psf,Nx,Ny)));
    %normalise
    Fpsf = Fpsf./max(Fpsf,[],'all');
    %Fourier transform ISM image
    Fimg = fftshift(fft2(data,Nx,Ny));
    fa = abs(Fimg); % modulus
    fp = angle(Fimg); % Phase angle
    %calculate reweighting factor
    Ffull = 1./(Fpsf + eps);
    Ffull(Fpsf < eps) = 0;
    %only reweight Fourier amplitude not phase
    W_ISM_img = abs(ifft2(Ffull.*fa.*exp(1i.*fp)));
    W_ISM_img = W_ISM_img(1:nx, 1:ny);
end