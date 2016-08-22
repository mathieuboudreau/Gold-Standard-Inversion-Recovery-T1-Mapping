% CustomFormat.m
%
% Sets axes properties
%
% written by J. Barral, M. Etezadi-Amoli, E. Gudmundson, and N. Stikov, 2009
%  (c) Board of Trustees, Leland Stanford Junior University 

axs = get(gcf, 'CurrentAxes');
myfigfont = 'Helvetica';

set( get(axs, 'Title'), 'FontName', myfigfont, 'FontSize', 24);
set( get(axs, 'Xlabel'), 'FontName', myfigfont, 'FontSize', 24);
set( get(axs, 'Ylabel'), 'FontName', myfigfont, 'FontSize', 24);
set( axs, 'FontName', myfigfont, 'FontSize', 22);
set( axs, 'Units', 'inches');