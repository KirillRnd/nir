function OM = OMfromAN(AN, k)
%OMfromAN вычисляет начальную фазу оптимального перелёта задачи встречи
% из угловой дальности
% OM (фаза) в витках, AN (угловая дальность) в витках
% k - соотношение радиусов

A=0.5433279;
C1 = A*(exp(k)+1)*exp(-k)*log(k);
C0 = -A*log(k)*exp(-k)*log(k);
OM = C1*AN+C0;
end