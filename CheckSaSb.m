rad = 1/32;
phi = 11/8;

%Первое задание переменных. Оно ломает мою программу
a_check = phi*pi - rad*pi;
b_check = phi*pi + rad*pi;

%Другое задание переменных. Тут всё работает как надо
a = (phi - rad)*pi;
b = (phi + rad)*pi;

%Вот тут возникает проблема
a_check-a
b_check-b