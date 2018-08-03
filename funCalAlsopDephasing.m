function sig=funCalAlsopDephasing(sig0, theta)
sig=sig0;
tmp=exp(theta*1i); a=real(tmp); b=imag(tmp);
sig(1)=sig0(1)*a-sig0(2)*b;
sig(2)=sig0(1)*b+sig0(2)*a;
end