
%% Buck_real %%

%close all
clear all

%parámetros de entrada:

n_bits_restric=12;  % Número de bits en las entradas y salidas
Admit_Rout=0.4;     % Admitancia de la Resistencia en la carga
Vg_value=12;        % Valor de Vg
dt=20e-9 ;          % diferencial de tiempo
L=0.000022 ;        % Inductuancia de la bobina en Henrios
C=0.000220;         % Capacitancia del condensador en Faradios
Vout_final= 5;      % Vout deseada        


%% inicializo variables: 
N_ciclos_M=2000;
N_ciclos_on_1=1;
N_ciclos_off_2=250;

N_ciclos_on_2=round(Vout_final*N_ciclos_off_2/Vg_value);
N_ciclos_off_1=N_ciclos_on_2+1;

i=1;
j=1;
il=0;
VoutAux=0;



%definimos Vg
% ------------------------------------------------
% Vg  Optimizado

Vg = Vg_value;  


%inicializamos Vout
Vout = 0 ; 
%inicializamos Vout_fb
Vout_fb = 0 ;

%dtl %
dtL=dt/L;

%dtC %
dtC=dt/C;


%% Fin de inicialización, inicio del bucle

for j=1:N_ciclos_M  % ciclos completos de reloj
for i=N_ciclos_on_1:N_ciclos_on_2  %mosfet on (parte de ciclo on)


% Vl Voltaje en la bobina dependiendo del Switch
Vl=(Vg-Vout_fb);
svLdtL=(Vl*dtL);
il=(il+svLdtL);
Iin=(il);  

%% Parte de abajo

ir=(Admit_Rout*Vout);
iL_fb=(il); 
iC=(iL_fb-ir);
siCdtC=(iC*dtC);
VoutAux=(VoutAux+siCdtC);

%%
Vout =(VoutAux);  
Vout_fb =(VoutAux); 

%% construir señales de salida.
VoutAux_senal(i,j)=VoutAux;
vl_senal(i,j)=Vl;
svLdtL_senal(i,j)=svLdtL;  
il_senal(i,j)=il;
iC_senal(i,j)=iC;
ir_senal(i,j)=ir;
siCdtC_senal(i,j)=siCdtC;
Vg_senal(i,j)=Vg;
Iin_senal(i,j)=Iin;

end

for i=N_ciclos_off_1:N_ciclos_off_2  %mosfet off (parte de ciclo off)

%Interruptor abierto

if il>0

    Vl=(-Vout_fb); 
    svLdtL=(Vl*dtL);
    il=(il+svLdtL);
    %

    Iin=0; % Diferente. igual que Iin pero valiendo 0
    %
    iL_fb=(il);   
    ir=(Admit_Rout*Vout);

    %
    iC=(iL_fb-ir); %No cambia
    
    %igual
    siCdtC=(iC*dtC);
    VoutAux=(VoutAux+siCdtC);
    Vout =(VoutAux);  
    Vout_fb =(VoutAux);


    
%% construir señales de salida.
VoutAux_senal(i,j)=VoutAux;
vl_senal(i,j)=Vl;
svLdtL_senal(i,j)=svLdtL; 
il_senal(i,j)=il;
iC_senal(i,j)=iC;
ir_senal(i,j)=ir;
siCdtC_senal(i,j)=siCdtC;
Vg_senal(i,j)=Vg;
Iin_senal(i,j)=Iin;

else 
    Vl=(0); %diferente
    svLdtL=(Vl*dtL);
    il=(il+svLdtL);
    Iin=(0); 
    iL_fb=(il);  
    ir=(Admit_Rout*Vout);

    iC=(-ir); % Detectada errata con respecto res_buck_completo_v0

    %igual
    siCdtC=(iC*dtC);
    VoutAux=(VoutAux+siCdtC);
    Vout =(VoutAux);  
    Vout_fb =(VoutAux);
    

%% construir señales de salida.
VoutAux_senal(i,j)=VoutAux;
vl_senal(i,j)=Vl;
svLdtL_senal(i,j)=svLdtL; 
il_senal(i,j)=il;
iC_senal(i,j)=iC;
ir_senal(i,j)=ir;
siCdtC_senal(i,j)=siCdtC;
Vg_senal(i,j)=Vg;
Iin_senal(i,j)=Iin;

end

% bucle para X ciclos de reloj con mosfet on y X ciclos con mosfet off.

end

end

%% Parte 2. Generación de señales real %%

VoutAux_real=reshape(VoutAux_senal,[],1); % convert matrix to column vector
Vl_real=reshape(vl_senal,[],1); % convert matrix to column vector
svLdtL_real=reshape(svLdtL_senal,[],1); % convert matrix to column vector
il_real=reshape(il_senal,[],1); % convert matrix to column vector
iC_real=reshape(iC_senal,[],1); % convert matrix to column vector
ir_real=reshape(ir_senal,[],1); % convert matrix to column vector
siCdtC_real=reshape(siCdtC_senal,[],1); % convert matrix to column vector
Vg_real=reshape(Vg_senal,[],1); % convert matrix to column vector
Iin_real=reshape(Iin_senal,[],1); % convert matrix to column vector


figure;
plot(VoutAux_real(:,1),'DisplayName','vout_senal(:,1)')


 %% Parte 3. Obtención automática de los valores en regimen permanente y transitório %%


 %ceil redondea hacia mas infinito, Floor hacia - infinito.

% Vg

Vg_max_trans=max_lejano_cero_trans(Vg_real);     % X
Vg_X= ceil(log2(Vg_max_trans))+1  ;

Vg_min_perm=mas_cercano_cero_perm(Vg_real,1000); % Y
Vg_Y= -(floor(log2(abs(Vg_min_perm))));

 
% vL

vL_max_trans=max_lejano_cero_trans(Vl_real)  ;   % X
vL_X= ceil(log2(vL_max_trans))+1  ; 

vL_min_perm=mas_cercano_cero_perm(Vl_real,1000); % Y
vL_Y= -(floor(log2(abs(vL_min_perm))));

    %vL_max_perm=mas_lejano_cero_perm(Vl_real,1000) % para los valores alternos
    %en permanente. No se utiliza por que no cruza por 0.

% svLdtL
%obtener max de svLdtL en transitorio.
 %svLdtL_max_trans= max(svLdtL_real);  % ----  
  svLdtL_max_trans=max_lejano_cero_trans(svLdtL_real);
  svLdtL_X= ceil(log2(svLdtL_max_trans))+1  ; 

 %obtener min de svLdtL en perm.
 % svLdtL_min_perm = min(svLdtL_real(length(svLdtL_real)-1000:length(svLdtL_real))); %HACER EN VALOR ABSOLUTO.
  svLdtL_min_perm=mas_cercano_cero_perm(svLdtL_real,1000);
  %svLdtL_max_perm=mas_lejano_cero_perm(svLdtL_real,1000) % para los valores alternos en permanente.
  svLdtL_Y= -(floor(log2(abs(svLdtL_min_perm))));


%il
il_max_trans=max_lejano_cero_trans(il_real)  ;   % X
il_X= ceil(log2(il_max_trans)) +1 ; 

il_min_perm=mas_cercano_cero_perm(il_real,1000); % Y
il_Y= -(floor(log2(abs(il_min_perm))));  %auqnue es triangular no cruza cero, se coge solo el valor mas pequeño.


%Iin
Iin_max_trans=max_lejano_cero_trans(Iin_real)  ;   % X
Iin_X= ceil(log2(Iin_max_trans))+1  ; 

%obtener max de Iin en perm.
Iin_max_perm=mas_lejano_cero_perm(Iin_real,1000);

Iin_min_perm=mas_cercano_cero_perm(Iin_real,1000); % Y

Iin_min_no_null=(5*(Iin_max_perm-Iin_min_perm)/100)/2; %al tener su minimo en cero he utilizado la regla cuando pasa por cero.

Iin_Y= -(floor(log2(abs(Iin_min_no_null)))); 


%iC
%obtener max de iC en transitorio.
iC_max_trans=max_lejano_cero_trans(iC_real)  ;   % X
iC_X= ceil(log2(iC_max_trans))+1  ; 

    % Como pasa por cero.
    %obtener max de iC en perm.
    iC_max_perm=mas_lejano_cero_perm(iC_real,1000);
    
    %obtener min de iC en perm.
    iC_min_perm=mas_cercano_cero_perm(iC_real,1000); % Y
    
    
    iC_min_no_null=(5*(iC_max_perm-iC_min_perm)/100)/2;

iC_Y= -(floor(log2(abs(iC_min_no_null))));


%ir
%obtener max de ir en transitorio.

ir_max_trans=max_lejano_cero_trans(ir_real)  ;   % X
ir_X= ceil(log2(ir_max_trans))+1  ; 

%obtener min de ir en permanente.
ir_min_perm=mas_cercano_cero_perm(ir_real,1000); % Y

ir_Y= -(floor(log2(abs(ir_min_perm))));

% ir_max_trans= max(ir_real);  % ----  

 %obtener max de siCdtC en transitorio.

siCdtC_max_trans=max_lejano_cero_trans(siCdtC_real)  ;   % X
siCdtC_X= ceil(log2(siCdtC_max_trans))+1  ; 

    % Como pasa por cero.
    %obtener max de siCdtC en perm.
    siCdtC_max_perm=mas_lejano_cero_perm(siCdtC_real,1000);
    
    %obtener min de siCdtC en perm.
    siCdtC_min_perm=mas_cercano_cero_perm(siCdtC_real,1000); % Y
    
    
    siCdtC_min_no_null=(5*(siCdtC_max_perm-siCdtC_min_perm)/100)/2;

siCdtC_Y= -(floor(log2(abs(siCdtC_min_no_null))));


 %ir
%obtener max de ir en transitorio.

VoutAux_max_trans=max_lejano_cero_trans(VoutAux_real)  ;   % X
VoutAux_X= ceil(log2(VoutAux_max_trans))+1  ; 

%obtener min de ir en permanente.
VoutAux_min_perm=mas_cercano_cero_perm(VoutAux_real,1000); % Y

VoutAux_Y= -(floor(log2(abs(VoutAux_min_perm))));

%crear señales dependientes.

% il_fb

il_fb_X=il_X;
il_fb_Y=il_Y;

% Vout_fb

Vout_fb_X=VoutAux_X;
Vout_fb_Y=VoutAux_Y;

% Vout

Vout_X=VoutAux_X;
Vout_Y=VoutAux_Y;

% Constantes
dtL_X= ceil(log2(dtL))  ; %ceil redondea hacia mas infinito, Floor hacia - infinito. No le pongo + 1 por se cte
dtL_Y= -(floor(log2(abs(dtL))));

dtC_X= ceil(log2(dtC))  ; %ceil redondea hacia mas infinito, Floor hacia - infinito. No le pongo + 1 por se cte
dtC_Y= -(floor(log2(abs(dtC))));



%% Agrupar
%% Accumulative

%igualo las Y
%n1v
siCdtC_Y;
VoutAux_Y;
if siCdtC_Y>VoutAux_Y
    VoutAux_Y=siCdtC_Y;
else
    siCdtC_Y=VoutAux_Y;
end

Y_Gr_n1v=siCdtC_Y;


%n1i
svLdtL_Y;
il_Y;
if svLdtL_Y>il_Y
    il_Y=svLdtL_Y;
else
    svLdtL_Y=il_Y;
end
Y_Gr_n1i=svLdtL_Y;





%% Non-accumulative
%n2v
Vg_Y;
Vout_fb_Y;
Vout_Y;
vL_Y;




G_n2v=[Vg_Y,Vout_fb_Y,Vout_Y,vL_Y];
Vg_Y=max(G_n2v);
Vout_fb_Y=max(G_n2v);
Vout_Y=max(G_n2v);
G_n2v;

Y_Gr_n2v=vL_Y;

%n2i

ir_Y;
il_fb_Y;
iC_Y;
Iin_Y;

G_n2i=[ir_Y,il_fb_Y,iC_Y,Iin_Y];
ir_Y=max(G_n2i);
il_fb_Y=max(G_n2i);
iC_Y=max(G_n2i);
Iin_Y=max(G_n2i);


Y_Gr_n2i=Iin_Y;



%% Constants

dtL_Y;
Y_Gr_n3L=dtL_Y;

dtC_Y;
Y_Gr_n3C=dtC_Y;

%% agrupar en 3
Vout_X;
il_X;


%comparar wl
%n1v y n1i
VoutAux_wl=VoutAux_X+VoutAux_Y+1; 
il_wl=il_X + il_Y +1;

if VoutAux_wl>il_wl
    n1i=VoutAux_wl-il_wl;
    il_Y=il_Y+n1i;
    svLdtL_Y=svLdtL_Y+n1i;
    n1v=0;
else
    n1v=il_wl-VoutAux_wl;
    %meter resto de resultados.
end

n1i;


%% En este grupo está la restricción de los bits. La entrada es Vg tiene que tener nº bits: (12)

%entrada Vg
n2v=n_bits_restric - Vg_X - Vg_Y;

Vg_Y_v2=Vg_Y+n2v;  %cambio para no sobrescribir
vL_Y_v2 =vL_Y+n2v;
Vout_fb_Y_v2=Vout_fb_Y+n2v;
Vout_Y_v2=Vout_Y+n2v;

%entrada Ir

n2i=n_bits_restric - ir_X - ir_Y;
ir_Y_v2=ir_Y+n2i;
il_fb_Y_v2=il_fb_Y+n2i;
iC_Y_v2=iC_Y+n2i;
Iin_Y_v2=Iin_Y+n2i;

%salida Iin. Ajuste de parte decimal de Iin para los 12 bits
Iin_X;
Iin_Y;

Iin_Y_v3 = n_bits_restric - Iin_X;

% salida Vout
Vout_Y_v3 = n_bits_restric - Vout_X;

%%  Comprobar wl de entradas, las del grupo n2v y n2i tener el mismo wl.
Vg_sfi=sfi(Vg,Vg_X+Vg_Y_v2,Vg_Y_v2);
Ir_sfi=sfi(ir,ir_X+ir_Y_v2,ir_Y_v2);

%salidas
Vout_sfi=sfi(Vout,Vout_X+Vout_Y_v3,Vout_Y_v3);
Iin_sfi=sfi(Iin,Iin_X+Iin_Y_v3,Iin_Y_v3);

%% Señales acumulativas. Igualar al peor caso. (realizado arriba)
VoutAux_sfi=sfi(VoutAux,VoutAux_X+VoutAux_Y+1,VoutAux_Y );
il_sfi=sfi(il,il_X+il_Y+1,il_Y );


%% añadir el n mas restrictivo a las ctes.
n1v;
n1i;
n2v;
n2i;

n_ctes=[n1v n1i n2v n2i];
n3L=max(n_ctes);
n3C=max(n_ctes);

dtL_sfi=sfi(dtL,dtL_X+dtL_Y+n3L,dtL_Y+n3L);  %sin bit para signo
dtC_sfi=sfi(dtC,dtC_X+dtC_Y+n3L,dtC_Y+n3C);  %sin bit para signo

dtL_Y_v2=dtL_Y+n3L;
dtC_Y_v2=dtC_Y+n3C;


%% mostrar resultados. Imprimir un fichero "number_bits.txt" con los resultados

salida.Non_accumulative.Vg_X=Vg_X;
salida.Non_accumulative.Vg_Y=Vg_Y_v2;

salida.Non_accumulative.vL_X=vL_X;
salida.Non_accumulative.vL_Y=vL_Y_v2;

salida.Accumulative.svLdtL_X=svLdtL_X;
salida.Accumulative.svLdtL_Y=svLdtL_Y;

salida.Accumulative.il_X=il_X;
salida.Accumulative.il_Y=il_Y;

salida.Non_accumulative.Iin_X=Iin_X;
salida.Non_accumulative.Iin_Y=Iin_Y_v3;

salida.Non_accumulative.il_fb_X=il_fb_X;
salida.Non_accumulative.il_fb_Y=il_fb_Y_v2;

salida.Non_accumulative.ir_X=ir_X;
salida.Non_accumulative.ir_Y=ir_Y_v2;

salida.Non_accumulative.iC_X=iC_X;
salida.Non_accumulative.iC_Y=iC_Y_v2;

salida.Accumulative.siCdtC_X=siCdtC_X;
salida.Accumulative.siCdtC_Y=siCdtC_Y;

salida.Accumulative.VoutAux_X=VoutAux_X;
salida.Accumulative.VoutAux_Y=VoutAux_Y;

salida.Non_accumulative.Vout_X=Vout_X;
salida.Non_accumulative.Vout_Y=Vout_Y_v3;

salida.Non_accumulative.Vout_fb_X=Vout_fb_X;
salida.Non_accumulative.Vout_fb_Y=Vout_fb_Y_v2;

salida.Constants.dtL_X=dtL_X;
salida.Constants.dtL_Y=dtL_Y_v2;

salida.Constants.dtC_X=dtC_X;
salida.Constants.dtC_Y=dtC_Y_v2;



text = jsonencode(salida,PrettyPrint=true);

fileID = fopen('number_bits.txt','w');
fprintf(fileID,'%s \n',text);
fclose(fileID);

%% funciones
function [resul]= mas_cercano_cero_perm (senal,n)
    a= max (senal(length(senal)-n:length(senal)));
    b= min(senal(length(senal)-n:length(senal)));
    aa=abs(a);
    bb=abs(b);
        if aa<bb
            resul=a;
        else
            resul=b;
        end
end


function [resul]= mas_lejano_cero_perm (senal,n)
    a= max (senal(length(senal)-n:length(senal)));
    b= min(senal(length(senal)-n:length(senal)));
    aa=abs(a);
    bb=abs(b);
        if aa>bb
            resul=a;
        else
            resul=b;
        end
end


function [resul]= max_lejano_cero_trans (senal)
    a= max (senal);
    b= min(senal);
    aa=abs(a);
    bb=abs(b);
        if aa>bb
            resul=a;
        else
            resul=b;
        end
end



