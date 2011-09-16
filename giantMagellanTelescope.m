classdef giantMagellanTelescope < telescopeAbstract
    
    properties (Constant)
        
        % segment diameter
        segmentD = 8.365;
        % radius of the ring where peripheric segments are located
        petalRingRadius = 8.658;
        % secondary mirror diameter
        centralObscurationD = 3.2;
        % center hole diameter 
        centerHoleD = 0;%0.55;
    end
    
    properties
        
        % # of segment
        nSegment = 7;
        % location of segment centers
        segmentCoordinate;
        % the gmt segment modeled as Zernike polynomials
        segment;
        tag = 'Giant Magellan Telescope';
        % the gmt segmented pupil
        pupil;
        
    end
    
    properties (Dependent)% , SetAccess = private)
        
        % diameter of a full telescope with the same area than gmt
        areaD;
        % pupil rotation angle        
        angle;
    end
    
    properties (Access=private)
        p_angle=0;
    end
        
    methods
        
        %% Constructor
        function obj = giantMagellanTelescope(varargin)
            
            p = inputParser;
            p.addParamValue('fieldOfViewInArcsec', [], @isnumeric);
            p.addParamValue('fieldOfViewInArcmin', [], @isnumeric);
            p.addParamValue('resolution', [], @isnumeric);
            p.addParamValue('samplingTime', [], @isnumeric);
            p.parse(varargin{:});

            % span diameter
            bigD = giantMagellanTelescope.petalRingRadius*2+giantMagellanTelescope.segmentD;
            obj = obj@telescopeAbstract(bigD,...
                'obstructionRatio',giantMagellanTelescope.centralObscurationD/giantMagellanTelescope.segmentD,...
                'fieldOfViewInArcsec',p.Results.fieldOfViewInArcsec,...
                'fieldOfViewInArcmin',p.Results.fieldOfViewInArcmin,...
                'samplingTime',p.Results.samplingTime,...
                'resolution',p.Results.resolution);
            setSegments(obj);
            display(obj);
        end
        
        %% Destructor        
        function delete(obj)
            cellfun(@(x)delete(x),obj.segment);
            checkOut(obj.log,obj)
        end
        
        %% Get telescope surface
        function out = area(obj)
            out = pi.*(7*giantMagellanTelescope.segmentD^2-6*giantMagellanTelescope.centerHoleD^2-giantMagellanTelescope.centralObscurationD^2)/4;
        end
        
        %% Get telescope diameter from full equivalente area 
        function out = get.areaD(obj)
            out = 2*sqrt(obj.area/pi);
        end
        
        %% Get/Set rotation angle
        function out = get.angle(obj)
            out = obj.p_angle;
        end
        function set.angle(obj,val)
            obj.p_angle = val;
            setSegments(obj);
        end
        
        %% Get and Set the pupil
%         function pupil = get.pupil(obj)
%             pupil = 0;
%             for k=1:obj.nSegment
%                 pupil = pupil + obj.segment{k}.pupil.*exp(1i.*obj.segment{k}.c);
%             end
%         end
        
        function display(obj)
            %% DISPLAY Display object information
            %
            % display(obj) prints information about the atmosphere+telescope object
            
%             display(obj.atm)
            fprintf('___ %s ___\n',obj.tag)
            fprintf(' %4.2fm diameter with a %4.2f%% central obstruction',...
                obj.D,obj.obstructionRatio*100)
            fprintf(' with %5.2fm^2 of light collecting area (%4.2fm diameter full single pupil);\n',obj.area,obj.areaD)
            if obj.fieldOfView~=0
                fprintf(' the field-of-view is %4.2farcmin;',...
                    obj.fieldOfView*constants.radian2arcmin)
            end
            if ~isempty(obj.resolution)
                fprintf(' the pupil is sampled with %dX%d pixels',...
                    obj.resolution,obj.resolution)
            end
            if obj.fieldOfView~=0 || ~isempty(obj.resolution)
                fprintf('\n')
            end
            fprintf('----------------------------------------------------\n')
        end
        
        function out = index(obj,z)
            %% INDEX
            
            out = false(size(z));
            for k=1:7
            out = out | ...
                ( abs(z -obj.segmentCoordinate(k)) < 0.5*obj.segmentD & ...
                abs(z -obj.segmentCoordinate(k)) > 0.5*obj.centerHoleD );
            end
            out = out & abs(z)>0.5*obj.centralObscurationD;
        end
        
        function out = otf(obj, zRho)
            %% OTF GMT optical transfert function
            %
            % out = otf(obj, x + iy) Computes the GMT optical transfert function
            
            out = diffPupAutoCorr(obj.segmentD,obj.centralObscurationD,abs(zRho));
                
            for k1 = 2:obj.nSegment
                r = zRho - obj.segmentCoordinate(k1);
                out = out + 2*diffPupCrossCorr(obj.segmentD,obj.centralObscurationD,obj.centerHoleD,abs(r));
                for k2 = 2:obj.nSegment
                     r = zRho - obj.segmentCoordinate(k1) - obj.segmentCoordinate(k2);
                     out = out + diffPupAutoCorr(obj.segmentD,obj.centerHoleD,abs(r));
                end
            end
            
            out = out/obj.area;
                        
%             if isa(obj.opticalAberration,'atmosphere')
%                 out = out.*phaseStats.otf(r,obj.opticalAberration);
%             end

            function out11 = diffPupAutoCorr(D,Dx,r)
                out11 = pupAutoCorr(D,r) + pupAutoCorr(Dx,r) - 2*pupCrossCorr(D/2,Dx/2,r);
            end

            function out21 = diffPupCrossCorr(D,Do,Dx,r)
                out21 = pupAutoCorr(D,r) - pupCrossCorr(D/2,Dx/2,r) - pupCrossCorr(D/2,Do/2,r) + pupCrossCorr(Do/2,Dx/2,r);
            end

            function out1 = pupAutoCorr(D,r)
                
                m_index       = r <= D;
                red         = r(m_index)./D;
                out1        = zeros(size(r));
                out1(m_index) = D.*D.*(acos(red)-red.*sqrt((1-red.*red)))./2;
                
            end
            
            function out2 = pupCrossCorr(R1,R2,r)
                
                out2 = zeros(size(r));
                
                m_index       = r <= abs(R1-R2);
                out2(m_index) = pi*min([R1,R2]).^2;
                
                m_index       = (r > abs(R1-R2)) & (r < (R1+R2));
                rho         = r(m_index);
                red         = (R1*R1-R2*R2+rho.*rho)./(2.*rho)/(R1);
                out2(m_index) = out2(m_index) + R1.*R1.*(acos(red)-red.*sqrt((1-red.*red)));
                red         = (R2*R2-R1*R1+rho.*rho)./(2.*rho)/(R2);
                out2(m_index) = out2(m_index) + R2.*R2.*(acos(red)-red.*sqrt((1-red.*red)));
                
            end

        end
        function out = psf(obj,fr,fo)  
            %% PSF GMT point spread function
            %
            % out = psf(obj, f) computes the GMT point spread function
            
%             if isa(obj.opticalAberration,'atmosphere')
%                 fun = @(u) 2.*pi.*quadgk(@(v) psfHankelIntegrandNested(v,u),0,obj.D);
%                 out = zeros(size(fr));
%                 parfor k = 1:numel(fr)
%                     out(k) = fun(f(k));
%                 end
% %                 out = arrayfun( fun, f);
%             else
                segmentArea    = pi.*giantMagellanTelescope.segmentD.^2./4;
                centerHoleArea = pi.*giantMagellanTelescope.centerHoleD.^2./4;
                out   = ones(size(fr)).*(segmentArea - centerHoleArea);
                index = fr~=0;
                u = pi.*giantMagellanTelescope.segmentD.*fr(index);
                v = pi.*giantMagellanTelescope.centerHoleD.*fr(index);
                out(index) = segmentArea.*2.*besselj(1,u)./u - centerHoleArea.*2.*besselj(1,v)./v;
                red = 0;
                for k=2:obj.nSegment
                    red = red + exp(-2.*1i.*pi.*fr.*abs(obj.segmentCoordinate(k)).*cos(fo-angle(obj.segmentCoordinate(k))) + 1i*obj.segment{k}.c);
                end
                out = out.*red;
                
%                 Ds = 1.5;
%                 C = obj.petalRingRadius;
%                 xL = C/2;
%                 yL = sqrt(3)*C/6;
%                 [oL,rL] = cart2pol(xL,yL);
%                 red = 0;
%                 for k = 1:6
%                     red = red + exp(-2.*1i.*pi.*fr(index).*rL.*cos(fo(index)-oL-(k-1)*pi/3));
%                 end
%                 u = pi.*Ds.*fr(index);
%                 surface = pi*Ds^2/4;
%                 out(index) = out(index) + red.*surface.*2.*besselj(1,u)./u;
% %                 out(~index) = out(~index) + red(~index).*surface;
%                 
%                 C = 2*obj.petalRingRadius;
%                 xL = C/2;
%                 yL = sqrt(3)*C/6;
%                 [oL,rL] = cart2pol(xL,yL);
%                 red = 0;
%                 for k = 1:6
%                     red = red + exp(-2.*1i.*pi.*fr(index).*rL.*cos(fo(index)-oL-(k-1)*pi/3));
%                 end
%                 u = pi.*Ds.*fr(index);
%                 surface = pi*Ds^2/4;
%                 out(index) = out(index) + red.*surface.*2.*besselj(1,u)./u;
% %                 out(~index) = out(~index) + red(~index).*surface;
                
                w = pi.*giantMagellanTelescope.centralObscurationD.*fr(index);
                centralObscurationArea = pi.*giantMagellanTelescope.centralObscurationD.^2./4;
                out(index) = out(index) + segmentArea.*exp(1i*obj.segment{1}.c).*2.*besselj(1,u)./u - centralObscurationArea.*2.*besselj(1,w)./w;
                out(~index) = out(~index) + segmentArea - centralObscurationArea;
                out = abs(out).^2./obj.area;
%             end
%             function y = psfHankelIntegrandNested(x,freq)
%                 y = x.*besselj(0,2.*pi.*x.*freq).*otf(obj,x);
%             end
        end
        
        function varargout = image(obj,resolution,pixelScaleInSpFreq)
            %% IMAGE 2D Point Spread Function
            %
            % psf = image(tel,resolution,pixelScaleInSpFreq) computes the
            % psf associated to the GMT pupil with the given resolution and
            % pixel scale
            %
            % image(tel,resolution,pixelScaleInSpFreq) displays the psf
            %
            % Example psf = image(gmt,128,arcsec(1e-3)/2.2e-2) computes the
            % psf with a resolution of 128x128 pixels and a pixel scale of
            % 1 mas at 2.2 micron
            
            n = resolution;
            u = pixelScaleInSpFreq*linspace(-1,1,n)*n/2;
            [fx,fy] = meshgrid(u);
            [fo,fr] = cart2pol(fx,fy);
            out = psf(obj,fr,fo)/psf(obj,0,0);
            
            if nargout>0
                varargout{1} = out;
            else
                imagesc(out)
                axis equal tight
            end
            
        end
        function out = fullWidthHalfMax(obj,fo)
            %% FULLWIDTHHALFMAX
            if nargin==1
                fo=0;
            end
            x0 = [0,2/obj.D];
            [out,~,exitflag] = fzero(@(x) psf(obj,abs(x)./2,fo) - psf(obj,0,fo)./2,x0,optimset('TolX',1e-9));
            if exitflag<0
                warning('cougar:telescope:fullWidthHalfMax',...
                    'No interval was found with a sign change, or a NaN or Inf function value was encountered during search for an interval containing a sign change, or a complex function value was encountered during the search for an interval containing a sign change.')
            end
            out = abs(out);
        end
        
        function obj = rmSegment(obj,segmentNumber)
            obj = giantMagellanTelescope(...
                'resolution',obj.resolution,...
                'fieldOfViewInArcsec',obj.fieldOfView*constants.radian2arcsec,...
                'samplingTime',obj.samplingTime);
            obj.nSegment = obj.nSegment - length(segmentNumber);
            obj.segment(segmentNumber)  = [];
            obj.segmentCoordinate(segmentNumber) = [];
        end

        function varargout = draw(obj,varargin)
            %% DRAW Draw pupil edges
            %
            % draw(gmt) draws the edges of the GMT segment
            %
            % draw(gmt,'PropertyName',propertyvalue,...) draws the edges of
            % the GMT segment using the values for the property
            % name/property value pairs specified for the lines
            %
            % h = draw(gmt,...) returns the handles of the graphic objects
            
            o = linspace(0,2*pi,101);
            co = cos(o);
            so = sin(o);
            segR = obj.segmentD/2; 
            obsR = obj.centralObscurationD/2; 
            handles = zeros(obj.nSegment+1,1);
            args = {'color','k','linewidth',2};
            if nargin>1
                args = varargin;
            end
            for kSegment = 1:obj.nSegment
                xc = real(obj.segmentCoordinate(kSegment));
                yc = imag(obj.segmentCoordinate(kSegment));
                handles(kSegment) = line(segR*co+xc,segR*so+yc,args{:});
            end
            handles(kSegment+1) = line(obsR*co,obsR*so,args{:});
            if nargout>0
                varargout{1} = handles;
            end
        end
        
    end
    
    methods (Access=private)
        function setSegments(obj)
            obj.segmentCoordinate = [ 0 , giantMagellanTelescope.petalRingRadius.*exp(1i.*(pi.*(0:5)/3 + obj.p_angle)) ];
            obstructionRatio = [ giantMagellanTelescope.centralObscurationD ones(1,6)*giantMagellanTelescope.centerHoleD]/giantMagellanTelescope.segmentD;
            for k=1:obj.nSegment
                [r,o] = utilities.cartAndPol(...
                    obj.resolution,...
                    obj.D/giantMagellanTelescope.segmentD,...
                    'offset',2*[real(obj.segmentCoordinate(k)),imag(obj.segmentCoordinate(k))]/giantMagellanTelescope.segmentD,...
                    'output','polar');
                obj.segment{k} = zernike(1,giantMagellanTelescope.segmentD,...
                    'obstructionRatio',obstructionRatio(k),...
                    'radius',r,'angle',o);
                obj.segment{k}.c = 0;
            end
            obj.pupil = 0;
            for k=1:obj.nSegment
                obj.pupil = obj.pupil + obj.segment{k}.pupil.*exp(1i.*obj.segment{k}.c);
            end
        end
    end
end