classdef Two_Photon_LIF < General_Operation_Point

    properties (Access=private)
    end
    
    methods (Access=public)
        function obj = Two_Photon_LIF(structinput,varargin)
            obj@General_Operation_Point(structinput,varargin);
            
        end
        
        function obj = calculate(obj)
            %% Keine Add-Ons notwendig! ff
            % Ausgangsdatei
            W1 = imread('Example.tif'); X1 = imresize(W1, 1);
            % figure
            % imshow(X1); title('1');
            % figure
            % imhist(X1)

            %% "IMAGE PROCESSING TOOLBOX" notwendig! ff
            % Kontrast verst�rken
            %precision=5;
            Y1 = histeq(X1); %Y1 = histeq(X1,precision);
            % figure
            % imshow(Y1); title('2');
            % figure
            % imhist(Y1)


            % Kontrast normen auf schwarz (scale max)
            Z1 = imadjust(Y1,stretchlim(Y1));
            figure
            imshow(Z1); title('Nozzle Koordinaten');
            %Nozzlekoordinaten input
            [x,y]=ginputRed(1);
            %Daten rechts des spiegelsymetrischen Sprays vernachl�ssigen
            Z2=Z1(:,1:floor(x));
            % figure
            % imhist(Z1)

            % Restkontrast normen auf wei� (scale min)
            A1 = Z2>0;
            % figure
            % imshow(A1); title('4');
            % figure
            % imhist(A1)

            %% "PARALLEL COMPUTING TOOLBOX" notwendig! ff
            % (Gewichteter und) Ungewichteter Schwerpunkt aller geschlossener Fl�chen
            data = regionprops(A1,Z2,{'Centroid','Area','EquivDiameter','MinorAxisLength','MajorAxisLength'}); %data = regionprops(A1,Z1,{'Centroid','WeightedCentroid'});
            % figure
            imshow(Z1); title('Ungewichteter Fl�chenschwerpunkte (blau), Nozzle (rot)'); %title('Gewichtet (rot) und Ungewichtete (blau) Fl�chenschwerpunkte'); 
            numObj = numel(data);
            hold on
            plot(x,y,'ro');
            for n = 1 : numObj
                % plot(data(n).WeightedCentroid(1), data(n).WeightedCentroid(2), 'ro')
                % Abstand zwischen Nozzle und Fl�chenschwerpunkt berrechnen und an struct �bergeben
                data(n).AbsoluteValue=sqrt(((data(n).Centroid(1)-x)^2)+(data(n).Centroid(2)-y)^2);
                %plot(data(n).Centroid(1), data(n).Centroid(2), 'b*')
            end
            hold off

            % Plot Pixelfl�che �ber Abstand Nozzle
            % figure
            % plot([data.AbsoluteValue],[data.Area],'o')
            figure
            %%%%%%HIER ALLGEMEINE REGRESSION MIT CHI-QUADRAT TEST!!
            f=polyfit([data.AbsoluteValue],[data.Area],2); x1=[data.AbsoluteValue];  y1=polyval(f,x1);
            plot(x1,y1,'b*'); title('Tr�pfenfl�che (-Volumen) �ber Nozzleanstand'); xlabel('Nozzleabstand [px]'); ylabel('Tr�pfchenfl�che [px^2]');
        end
    end
end

