function atm = gmtAtmosphere(ID,varargin)

p = inputParser;
p.addRequired('ID', @isnumeric);
p.addParamValue('L0',60,@isnumeric);
p.parse(ID, varargin{:});
ID = p.Results.ID;
L0 = p.Results.L0;

z  = [25 275 425 1250 4000 8000 13000];
fprintf(' *** GMT ATMOSPHERE \n')
switch ID
    case 1
        fprintf('Case 1: MEDIAN (Goodwin LCO Jan08, case 5)\n ***\n')
        r0 = 0.151;
        fr0 = [0.1257 0.0874 0.0666 0.3498 0.2273 0.0681 0.0751];
        vs = [5.6540 5.7964 5.8942 6.6370 13.2925 34.8250 29.4187];
        vd = [0.7798 8.2564 12.4752 32.5007 72.0987 93.1994 100.0507]*pi/180;
        tag = 'GMT-ATM1';
    case 2
        fprintf('Case 2: STRONG GROUND LAYER (Goodwin LCO Jan08, case 7)\n *** \n')
        r0 = 0.151;
        fr0 = [0.2039 0.1116 0.1706 0.2281 0.1694 0.0615 0.0549];
        vs = [9.4092 9.6384 9.7905 2.7028 8.2883 30.1824 13.0364];
        vd = [0.7798 8.2564 12.4752 32.5007 72.0987 93.1994 100.0507]*pi/180;
        tag = 'GMT-ATM2';
    case 3
        fprintf('Case 3: GOOD CONDITIONS (Goodwin LCO Jan08, case 1)\n *** \n')
        r0 = 0.19;
        fr0 = [0.1410 0.0381 0.0542 0.3403 0.2528 0.0917 0.0819];
        vs = [2.1953 2.2575 2.3029 2.7028 8.2883 30.1824 13.0364];
        vd = [0.7798 8.2564 12.4752 32.5007 72.0987 93.1994 100.0507]*pi/180;
        tag = 'GMT-ATM3';
    case 4
        fprintf('Case 4: POOR CONDITIONS (Goodwin LCO Jan08, case 9)\n *** \n')
        r0 = 0.12;
        fr0 = [0.1335 0.0730 0.1117 0.3311 0.2170 0.0613 0.0724];
        vs = [9.4092 9.6384 9.7905 10.8527 18.2550 39.1520 43.7936];
        vd = [0.7798 8.2564 12.4752 32.5007 72.0987 93.1994 100.0507]*pi/180;
        tag = 'GMT-ATM4';
    case 5
        fprintf('Case 5: THIN CLOUDS (Goodwin LCO Jan08, case 5)\n *** \n')
        r0 = 0.12;
        fr0 = [0.1335 0.0730 0.1117 0.3311 0.2170 0.0613 0.0724];
        vs = [9.4092 9.6384 9.7905 10.8527 18.2550 39.1520 43.7936];
        vd = [0.7798 8.2564 12.4752 32.5007 72.0987 93.1994 100.0507]*pi/180;
        tag = 'GMT-ATM5';
    otherwise
        error('Oups! only 5 cases, so far!')
end

atm = atmosphere(500e-9,r0,L0,...
    'altitude',z,...
    'fractionnalR0',fr0,...
    'windSpeed',vs,...
    'windDirection',vd);
atm.tag = tag;