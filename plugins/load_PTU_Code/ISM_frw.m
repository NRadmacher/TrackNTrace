function W_ISM_img = ISM_frw(data, psf)    
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
    psf = psf - min(psf,[],'all');
    %remove possible NaNs
    psf(isnan(psf)) = 0.0;
    %Fouier transform
    Fpsf = abs(fftshift(fft2(psf,Nx,Ny)));
    %normalise
    Fpsf = Fpsf./max(Fpsf,[],'all');

    %regulation parameter to prevent divergence
    eps = 0.1;
    %Fourier transform ISM image
    Fimg = fftshift(fft2(data,Nx,Ny));
    fa = abs(Fimg); % modulus
    fp = angle(Fimg); % Phase angle
    %calculate reweighting factor
    Ffull = 1./(Fpsf + eps);
    Ffull(Fpsf < eps) = 0;
    %only reweight Fourier amplitude not phase
    W_ISM_img = abs(ifftshift(ifft2(ifftshift(Ffull.*fa.*exp(1i.*fp)))));
    W_ISM_img = W_ISM_img(Nx/2:Nx/2+nx-1, Ny/2:Ny/2+ny-1);
end