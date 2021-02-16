function[u] = conservative(t,ak,phase,w)
% generate input signal
for k=1:length(phase)
    fun(k) = ak*cos(2*pi*w(k)*t + phase(k)); 
end
u = sum(fun);