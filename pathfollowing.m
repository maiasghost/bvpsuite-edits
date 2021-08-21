function [coeff,x1,valx1,x1tau,valx1tau,polynomials]=pathfollowing(bvpfile,zeichnung,x1,start,ersterschritt,schrittanzahl,bvpopt,ausgabe,pfad,startindex,speichername,pfaddatenhelp,gitter,feinesgitter,hitpoint)
%Soll bereits mit einer ausgerechneten L�sung aufgerufen werden
%Dieses Ergebnis (start) sollte den unbekannten Parameter p1 an der Stelle
%N(nm+summe(ordnung))+1 haben
%Die Gleichung, die den Parameter spezifiziert (z.B. p1=7) mu� unbedingt im
%Anschlu� an die Rand- bzw. Zusatzbedingungen stehen, da sie genau an
%dieser Stelle ersetzt wird. (Also die Stelle
%summe(ordnung)+1)

%bvpfile: Datei, die allgemeine Daten �ber das Problem enth�lt (im �blichen
%Format)
%zeichnung: bin�r, ob nach jeder Berechnung eine Graphik angefertigt wird
%x1: Gitter der Startl�sung
%start: Koeffizienten der Startl�sung
%ersterschritt: gr��e des ersten Schrittes bez�glich des Parameters, die
%weiteren Tangen haben dieselbe l�nge, aber i.a. andere Richtungen, was
%bewirkt, dass die n�chsten Schritte bez�glich des Parameters nicht mehr
%gleich weit sind
%schrittanzahl: maximale Anzahl von Durchl�ufen der Schleife
%bvpopt: options bez�glich Newton-Solver usw. im Standard-Format
%ausgabe: Standard Parameter, ob w�hrend der Berechnung Ausgaben get�tigt
%werden
%pfad: Pfad, wo gespeichert wird
%startindex: Gibt an, an welcher Stelle des Pfades die Berechnung
%fortgesetzt wird
%pfaddatenhelp: Geben Sie hier an, welche signifikanten Dinge der L�sung
%gesondert gespeichert werden sollen, um daraus dann einen Pfad graphisch
%darzustellen. Es handelt sich dabei um eine ? x 2 Matrix (? = Anzahl der
%darzustellenden Pfade). Am Beginn jeder Zeile steht, welche
%L�sungskomponente Sie darstellen wollen. 0 codiert die Parameter, die in
%der Gleichung vorkommen. Der 2. Teil gibt an, welche Stelle
%herausgegriffen werden soll.
%z.B. [3,0.4;0,2;1,0] speichert 3 Pfade, der erste stellt den Parameter
%gegen z3(0.4) dar, der 2. den Parameter gegen p2 und der 3. den Parameter
%gegen z1(0). Geben Sie den Spaltenvektor(!) [?;1] an, wenn Sie die
%Maximumnorm der Komponente ? w�nschen
%feinesgitter: boolean, gibt an, ob bei der Gittersteuerung die L�sung auch
%auf dem feinen Gitter ausgegeben werden soll

%Will man das Programm mit der Gittersteuerung verwenden, so m�ssen
%folgende Parameter in der Datei options.mat gespeichert sein:
%abstolgitter,reltolgitter,K,wiederholung

%hitpoint: A value for p1 which will be hit in the path - to be implemented
%in future versions


if nargin<15 hitpoint=[]; if nargin<14 feinesgitter=1; if nargin<13 gitter=0; if nargin<12 pfaddatenhelp=[]; if nargin<11 speichername='path.mat'; if nargin<10 startindex=1; if nargin<9 pfad=[]; if (nargin<8) ausgabe=1; if (nargin<7) bvpopt=[]; if (nargin<6) schrittanzahl=1; if (nargin<5) schrittweite=0.1; if (nargin<4) start=[];if(nargin<3) x1=[];if(nargin<2) zeichnung=false;end;end;end;end;end;end;end;end;end;end;end;end;end;end;
% global DF help1
%Nummer des Parameters, der verfolgt werden soll (im Moment nur f�r 1 implementiert, da equations.m angepa�t werden mu�):
parnum=1;
schrittweite=0; %schrittweite is just an internal variable ... the real stepsize is saved in ersterschritt
praediktor=[];

%Definieren aller bvpfile-spezifischen Parameter:

ordnung=feval(bvpfile,'ordnung');
n=length(ordnung);
anz_parameter=feval(bvpfile,'parameter');
%Berechnung des Tangentenvektors auf die L�sung
schritt=startindex;
par=feval(bvpfile,'parameter');

%tempor�rer Test!
%hitpoint=0.12;


x1tau=equations('x1tau',start,bvpfile,x1);
try
    valx1=equations('valx1',start,bvpfile,x1);
    valx1tau=equations('valx1tau',start,bvpfile,x1);
catch
    valx1=0;
    valx1tau=0;
end
if startindex>1
    load(strcat(pfad,speichername));
end
parametervalue(startindex)=start(end-anz_parameter+parnum);


if length(pfaddatenhelp(1,:))==2
    for i=1:length(pfaddatenhelp(:,1))
        if pfaddatenhelp(i,1)~=0
            hilf=equations('wert',start,bvpfile,x1,pfaddatenhelp(i,2));
            pathdata(i,startindex)=hilf(pfaddatenhelp(i,1));
        else
            pathdata(i,startindex)=start(end-anz_parameter+pfaddatenhelp(i,2));
        end
    end
else
    pathdata(1,startindex)=max(abs(valx1tau(pfaddatenhelp(1,1),:)));
end

%pathdata(startindex)=valx1tau(2,1);

speichername2=speichername;
if length(strfind(speichername,'.mat'))>0
    speichername2=speichername(1:strfind(speichername,'.mat')-1);
end

profilname=strcat(strcat(speichername2,num2str(startindex)),'.mat');
coeff=start;
save(strcat(pfad,profilname),'coeff','x1','valx1','x1tau','valx1tau');

x1last=[];
startlast=[];
praediktor(1,:)=start.';
passed_hitpoint=0; %Boolean variable to figure out if the hitpoint is passed or not
if length(hitpoint)>0 %Test, if the parametervalue is at the moment smaller or bigger than hitpoint
    if parametervalue(startindex)<hitpoint
        smaller_hitpoint=1;
    else
        smaller_hitpoint=0;
    end
end
halvestep=0; %If stepsize is too long, halvestep will be set to 1 and the stepsize will be halved;
ersterschritt_original=ersterschritt; %save original stepsize endless=0;
schrittlast=schritt;
while schritt<startindex+schrittanzahl
        if ~halvestep
            if ~passed_hitpoint
                schritt=schritt+1;
            else
                profilname=strcat(strcat(speichername2,num2str(schritt-1)),'.mat');
                lastfile=strcat(pfad,profilname);
                load (lastfile);
                start=coeff;
            end
        else
            halvestep=0;
            start=start_halve;
            praediktor=praediktor_halve;
            x1=x1_halve;
            schritt=schritt+1;
            if schritt==schrittlast
                 endless=endless+1;
            else
                endless=0;
            end;
            schrittlast=schritt;
            if schritt==startindex+1
                ersterschritt=ersterschritt/2;
            else
                schrittweite=schrittweite/2;
            end;
            if endless>100
                error('The maximum number of step halvings is reached!');
            end
            fprintf('The step is too big and would cause an error. The stepsize has been halved!\n');
        end
        try
        start_halve=start; %save all original variables for necessary halving;
        praediktor_halve=praediktor;
        x1_halve=x1;
        if ~passed_hitpoint
            [start,tangente,schrittweite]=tangente_berechnen(ordnung,ersterschritt,schritt,startindex,parametervalue,start,schrittweite,bvpfile,praediktor,x1);
            if schritt==startindex+1
                schrittweite_original=schrittweite;
            end
        else
            %When the hitpoint is passed at the first step, the stepsize
            %should not be recalculated
            [start,tangente,schrittweite]=tangente_berechnen(ordnung,ersterschritt,startindex+2,startindex,parametervalue,start,schrittweite,bvpfile,praediktor,x1);
        end
        tangentebak(schritt,1:length(tangente))=tangente;
        if ~passed_hitpoint
            if sign(tangentebak(schritt-1,:)*tangente)<0
                fprintf('Found turning point near %f!\n',parametervalue(schritt-1));
                schrittweite=-schrittweite;
            end
        end
        start=start+schrittweite*tangente;
        praediktor(1,1:length(start))=start.';
        praediktor(2,1:length(tangente))=tangente.';
        if gitter==0
        [start,x1,valx1,x1tau,valx1tau,polynomials]=run(bvpfile,zeichnung,x1,start,bvpopt,ausgabe,praediktor);
        end
        if gitter
        load options abstolgitter reltolgitter K wiederholung
        [start,x1,valx1,x1tau,valx1tau,polynomials]=run(bvpfile,zeichnung,x1,start,bvpopt,ausgabe,praediktor);
        [start_mesh,x1_mesh,valx1_mesh,x1tau_mesh,valx1tau_mesh,polynomials_mesh,valerror_mesh,toleriert_mesh]=meshadaptation(abstolgitter,reltolgitter,K,bvpfile,zeichnung,x1,start,bvpopt,ausgabe,wiederholung,feinesgitter,praediktor);
        end


    profilname=strcat(strcat(speichername2,num2str(schritt)),'.mat');
    coeff=start;
    if gitter==0
        save(strcat(pfad,profilname),'coeff','x1','valx1','x1tau','valx1tau');
    end
    if gitter
        coeff_bak=coeff;
        x1_bak=x1;
        valx1_bak=valx1;
        x1tau_bak=x1tau;
        valx1tau_bak=valx1tau;
        coeff=start_mesh;
        x1=x1_mesh;
        valx1=valx1_mesh;
        x1tau=x1tau_mesh;
        valx1tau=valx1tau_mesh;
        save(strcat(pfad,profilname),'coeff','x1','valx1','x1tau','valx1tau');
        coeff=coeff_bak;
        x1=x1_bak;
        valx1=valx1_bak;
        x1tau=x1tau_bak;
        valx1tau=valx1tau_bak;
    end
    if ~gitter
        parametervalue(schritt)=start(end-anz_parameter+parnum);
    else
        parametervalue(schritt)=start_mesh(end-anz_parameter+parnum);
    end

    if length(pfaddatenhelp(1,:))==2
        for i=1:length(pfaddatenhelp(:,1))
            if pfaddatenhelp(i,1)~=0
                if ~gitter
                    hilf=equations('wert',start,bvpfile,x1,pfaddatenhelp(i,2));
                else
                    hilf=equations('wert',start_mesh,bvpfile,x1_mesh,pfaddatenhelp(i,2));
                end
                pathdata(i,schritt)=hilf(pfaddatenhelp(i,1));
            else
                pathdata(i,schritt)=start(end-anz_parameter+pfaddatenhelp(i,2));
            end
        end
    else
        pathdata(1,schritt)=max(abs(valx1tau(pfaddatenhelp(1,1),:)));
    end
    if length(pfaddatenhelp)>0
        save(strcat(pfad,speichername),'parametervalue','pathdata');
    else
        save(strcat(pfad,speichername),'parametervalue');
    end
    fprintf('Parametervalue at index %i: %f\n',schritt,parametervalue(schritt));

        if length(hitpoint)>0
            if ~passed_hitpoint
                if abs(parametervalue(schritt)-hitpoint)<1e-5
                    fprintf('Hitpoint found at index %i\n',schritt);
                    passed_hitpoint=0;
                    schrittweite=sign(schrittweite)*abs(schrittweite_original);
                    freestep=1;
                    if smaller_hitpoint
                        smaller_hitpoint=0;
                    else
                        smaller_hitpoint=1;
                    end
                else
                    if smaller_hitpoint
                        if parametervalue(schritt)>hitpoint
                            passed_hitpoint=1;
                            fprintf('Hitpoint passed at index %i\n',schritt);
                            schrittweite=schrittweite/2;
                        end
                    else
                        if parametervalue(schritt)<hitpoint
                            passed_hitpoint=1;
                            fprintf('Hitpoint passed at index %i\n',schritt);
                            schrittweite=schrittweite/2;
                        end
                    end
                end
            else %If hitpoint has been passed
               if abs(parametervalue(schritt)-hitpoint)<1e-5
                    fprintf('Hitpoint found at index %i\n',schritt);
                    passed_hitpoint=0;
                    schrittweite=sign(schrittweite)*abs(schrittweite_original);
                    freestep=1;
                    if smaller_hitpoint
                        smaller_hitpoint=0;
                    else
                        smaller_hitpoint=1;
                    end
                else
                    if smaller_hitpoint
                        if parametervalue(schritt)<hitpoint
                            passed_hitpoint=0;
                            schrittanzahl=schrittanzahl+1;
                        else
                            passed_hitpoint=1;
                            fprintf('Hitpoint passed at index %i\n',schritt);
                            schrittweite=schrittweite/2;
                        end
                    else
                        if parametervalue(schritt)>hitpoint
                            passed_hitpoint=0;
                            schrittanzahl=schrittanzahl+1;
                        else
                            passed_hitpoint=1;
                            fprintf('Hitpoint passed at index %i\n',schritt);
                            schrittweite=schrittweite/2;
                        end
                    end
                end
            end
        end
        catch
            halvestep=1;
            schritt=schritt-1;
        end


end
fprintf('Calculation finished.\n');
coeff=start;


function [start,tangente,schrittweite]=tangente_berechnen(ordnung,ersterschritt,schritt,startindex,parametervalue,start,schrittweite,bvpfile,praediktor,x1);


    DF=equations('DF',start,bvpfile,x1);
    %Berechnen des Tangentenvektors--------------------
    %DF ist eine quadratische Matrix. Die Zeile, die die Ableitung der
    %Gleichung, wo p1 spezifiziert wurde, enth�lt, mu� "weggedacht" werden, um eine
    %Matrix mit einem um 1 verminderten Rang zu erhalten und daher einen
    %Tangentenvektor berechnen zu k�nnen

    %Da DF vollen Rang hat wird dies aber gleich ausgenutzt und das lineare
    %Gleichungssystem DF*tangente=[0,0,...,0,1,0,...,0] gel�st, dabei steht
    %der
    %Einser an der summe(ordnung)+1-ten Stelle. tangente ist also ein
    %Tangentenvektor zur um eine Zeile verminderten Matrix von DF
    help1=0;
    help1(1:sum(ordnung),1)=0;
    help1(sum(ordnung)+1,1)=abs(ersterschritt);
    if length(DF(:,1))>sum(ordnung)+1
        help1(sum(ordnung)+2:length(DF(:,1)),1)=0;
    end

    tangente=DF\help1;

    if schritt==startindex+1
        schrittweite=sign(ersterschritt)*sqrt(sum(tangente.*tangente));
    end
    tangente=tangente/sqrt(sum(tangente.*tangente));