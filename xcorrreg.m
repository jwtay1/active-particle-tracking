function pxShift = xcorrreg(refImg, movedImg)

%REGISTERIMG  Register two images using cross-correlation
%
%  I = xcorrreg(R, M) registers two images by calculating the
%  cross-correlation between them. R is the reference or stationary image,
%  and M is the moved image.
%
%  Note: This algorithm only works for translational shifts, and will not
%  work for rotational shifts or image resizing.

%Compute the cross-correlation of the two images
crossCorr = ifft2((fft2(refImg) .* conj(fft2(movedImg))));

%Find the location in pixels of the maximum correlation
[xMax, yMax] = find(crossCorr == max(crossCorr(:)));

%Compute the relative shift in pixels
Xoffset = fftshift(-size(refImg,1)/2:(size(refImg,1)/2 - 1));
Yoffset = fftshift(-size(refImg,2)/2:(size(refImg,2)/2 - 1));

pxShift = round([Xoffset(xMax), Yoffset(yMax)]);

end