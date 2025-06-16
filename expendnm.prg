closeall;
clearg NJ,TJ,indy,indx,indz,indiv,w,das,sco,hes,
       datn,datt,N,indyear,minyear,maxyear,years,w2,tmax,ji,
       zgy,h,justid,varx,varz,__title,mliter,staf;

format 10,4;

output off;

proc makeyears;
    local years,dt,j,sji,ji,yx,year;


    years = 0;

    open dt = ^dataset;
    j = 1;
    do until j>rows(datn);
        NJ = datn[j];
        TJ = datt[j];

        sji = 1;
        do until sji>NJ;
            year = readr(dt,TJ);
            if TJ<lagl+llev+1;
                goto jy;
            endif;
            year = year[.,indyear];
            years = years|year[1+lagl:rows(year)];
            years = years[uniqindx(years,1)];
jy:
            sji = sji + 1;
        endo;
        j = j+1;
    endo;
    dt = close(dt);
    retp(trimr(years,1,0));
endp;


proc makez(q,l1,l2);
    local j,z,dlag;

    q = q';
    j = 1;
    das = eye(tmax-lagl);

    dlag = l1-lagl;
    if dlag>0;
        j = dlag+1;
        das = eye(tmax-l1);
    endif;
    dlag = dlag*(dlag>0);
    z = q[maxc(1|lagl-l2+1+dlag):minc((lagl-l1+1+dlag)|cols(q))].*.das[.,1];
    j = 2;
    do until j>cols(das);
        z = z~(q[maxc(1|lagl-l2+j+dlag):minc((lagl-l1+j+dlag)|cols(q))]
                .*.das[.,j]);
        j = j+1;
    endo;
    if dlag>0;
        z = zeros(dlag,cols(z))|z;
    endif;
    retp(z);
endp;


proc makezi(yx,year,ji);
    local ziv,jziv,i,per,qi,perf,zi,jz,z;

    if nonseqz>0;
        ziv = lagn(yx[.,indiv[1]],lnseqz[1]);
        jziv = 2;
        do until jziv>rows(indiv);
            ziv = ziv~lagn(yx[.,indiv[jziv]],lnseqz[jziv]);
            jziv = jziv+1;
        endo;
    endif;

    i = 1;
    do until i>ji;
        if seqz>0;
            per = year[(i-1)*TJ+1]-minyear;
            qi = yx[(i-1)*TJ+1:i*TJ,indz];
            if per>0;
                qi = zeros(per,cols(qi))|qi;
            endif;
            perf = tmax-rows(qi);
            if perf>0;
                qi = qi|zeros(perf,cols(qi));
            endif;
            zi = makez(qi[.,1],lseqz1[1],lseqz2[1]);
            jz = 2;
            do until jz>rows(indz);
                zi = zi~makez(qi[.,jz],lseqz1[jz],lseqz2[jz]);
                jz = jz+1;
            endo;

            zi = trimr(zi,per,perf);
            if nonseqz>0;
                zi = zi~(ziv[(i-1)*TJ+lagl+1:i*TJ,.]);
            endif;
        elseif nonseqz>0;
            zi = ziv[(i-1)*TJ+lagl+1:i*TJ,.];
        endif;
        if timedum==1 or timez==1;
            if seqz>0 or nonseqz>0;
                zi = (year[(i-1)*TJ+lagl+1:i*TJ].==years')~zi;
            else;
                zi = (year[(i-1)*TJ+lagl+1:i*TJ].==years');
            endif;
            if model==1;
                zi = trimr(zi',1,0)';
            endif;
        elseif model==0;
            zi = ones(rows(zi),1)~zi;
        endif;
        if i==1;
            z = zi;
        else;
            z = z|zi;
        endif;
        i = i+1;
    endo;
    retp(z);
endp;


proc makezz;
    local zz,dt,j,sji,ji,yx,year,z;


    zz = 0;

    open dt = ^dataset;
    j = 1;
    do until j>rows(datn);
        NJ = datn[j];
        TJ = datt[j];

        sji = 0;
        ji = minc(nind|NJ);
        do until sji==NJ;
            yx = readr(dt,ji*TJ);
            if TJ<lagl+llev+1;
                goto jpp;
            endif;
            year = yx[.,indyear];
            z = makezi(yx,year,ji);
            zz = zz + z'z;
jpp:
            sji = sji + ji;
            ji = minc(nind|(NJ-sji));
        endo;
        j = j+1;
    endo;
    dt = close(dt);
    retp(zz);
endp;

proc makezuz(bb);
    local a,b,zz,dt,j,sji,ji,yx,mom,z,db,i,momi,zi;

    if lfmy>0;
        a = bb[1:lfmlag];
        b = bb[lfmlag+1:rows(bb)];
    else;
        a = 0;
        b = bb;
    endif;

    zz = 0;

    open dt = ^dataset;
    j = 1;
    do until j>rows(datn);
        NJ = datn[j];
        TJ = datt[j];
        sji = 0;
        ji = minc(nind|NJ);
        do until sji==NJ;
            yx = readr(dt,ji*TJ);
            if TJ<lagl+llev+1;
                goto jp;
            endif;
            {mom,z,db} = makemom(yx,ji,j,a,b);
            i = 1;
            do until i>ji;
                momi = mom[(i-1)*(TJ-lagl)+1:i*(TJ-lagl)];
                zi = z[(i-1)*(TJ-lagl)+1:i*(TJ-lagl),.];
                zz = zz + zi'momi*momi'zi;
                i = i+1;
            endo;
jp:
            sji = sji+ji;
            ji = minc(nind|(NJ-sji));
        endo;
        j = j+1;
    endo;
    dt = close(dt);
    retp(zz);
endp;



proc (3)=makemom(yx,ji,j,a,b);
    local y,x,z,jx,year,yeare,e,e1,y1,x1,
          ly,jy,ly1,mu,mu1,mu1mu,mom,da,db,
          my,mly,jly,mmu,mmux,jmx;

    y = yx[.,indy];
    x = lagn(yx[.,indx[1]],lx[1]);
    jx = 2;
    do until jx>rows(indx);
        x = x~lagn(yx[.,indx[jx]],lx[jx]);
        jx = jx+1;
    endo;

    year = yx[.,indyear];
    yeare = year;

    if lagl>0;
        e = ones(TJ,1);
        e[1:lagl] = zeros(lagl,1);
        e = ones(ji,1).*.e;

        if model==2;
            e1 = ones(TJ,1);
            if lagl>1;
                e1[1:lagl-1] = zeros(lagl-1,1);
            endif;
            e1[datt[j]] = 0;
            e1 = ones(ji,1).*.e1;
            y1 = selif(y,e1);
            x1 = selif(x,e1);
        endif;
        y = selif(y,e);
        x = selif(x,e);
        yeare = selif(year,e);
    endif;

    if lfmy>0;
        ly = lag1(yx[.,indy]);
        jy = 2;
        do until jy>lfmlag;
            ly = ly~lagn(yx[.,indy],jy);
            jy = jy+1;
        endo;
        if model==2;
            ly1 = selif(ly,e1);
        endif;
        ly = selif(ly,e);
    endif;


    if timedum==1;
        x = (yeare.==years')~x;
        if model==2;
            x1 = zeros(rows(x1),rows(years))~x1;
        elseif model==1;
            x = trimr(x',1,0)';
        else;
            x[.,1] = ones(rows(x),1);
        endif;
    elseif model==0;
        x = ones(rows(x),1)~x;
    endif;

    z = makezi(yx,year,ji);

    mu = exp(x*b);
    if model==2;
        mu1 = exp(x1*b);
        mu1mu = exp((x1-x)*b);
    endif;

    if lfmy>0;
        if model==2;
            if qdif==1;
                mom = ((y-ly*a).*mu1mu)-(y1-ly1*a);
                da = ly1-ly.*mu1mu;
                db = ((y-ly*a).*mu1mu).*(x1-x);
            else;
                mom = ((y-ly*a)./mu)-((y1-ly1*a)./mu1);
                da = (ly1./mu1)-(ly./mu);
                db = (((y1-ly1*a)./mu1).*x1)-(((y-ly*a)./mu).*x);
            endif;
        elseif model==1;
            my = meanc(reshape(y,ji,TJ-lagl)').*.ones(TJ-lagl,1);
            mly = meanc(reshape(ly[.,1],ji,TJ-lagl)');
            jly = 2;
            do until jly>cols(ly);
                mly = mly~(meanc(reshape(ly[.,jly],ji,TJ-lagl)'));
                jly = jly+1;
            endo;
            mly = mly.*.ones(TJ-lagl,1);
            mmu = meanc(reshape(mu,ji,TJ-lagl)').*.ones(TJ-lagl,1);
            mmux = meanc(reshape(mu.*x[.,1],ji,TJ-lagl)');
            jmx = 2;
            do until jmx>cols(x);
                mmux = mmux~meanc(reshape(mu.*x[.,jmx],ji,TJ-lagl)');
                jmx = jmx+1;
            endo;
            mmux = mmux.*.ones(TJ-lagl,1);
            mom = y-ly*a-mu.*(my-mly*a)./mmu;

            da = -ly + mu.*mly./mmu;
            db = -(mu.*(my-mly*a)./mmu).*x
                       + (mu.*(my-mly*a)./(mmu^2)).*mmux;
        else;
            if addit==1;
                mom = y-ly*a-mu;
                da = -ly;
                db = -mu.*x;
            else;
                mom = ((y-ly*a)-mu)./mu;
                da = -ly./mu;
                db = -((y-ly*a)./mu).*x;
            endif;
        endif;
        db = da~db;
    else;
        if model==2;
            if qdif==1;
                mom = (y.*mu1mu)-y1;
                db = y.*mu1mu.*(x1-x);
            else;
                mom = (y./mu)-(y1./mu1);
                db = ((y1./mu1).*x1)-((y./mu).*x);
            endif;
        elseif model==1;
            my = meanc(reshape(y,ji,TJ-lagl)').*.ones(TJ-lagl,1);
            mmu = meanc(reshape(mu,ji,TJ-lagl)').*.ones(TJ-lagl,1);
            mmux = meanc(reshape(mu.*x[.,1],ji,TJ-lagl)');
            jmx = 2;
            do until jmx>cols(x);
                mmux = mmux~meanc(reshape(mu.*x[.,jmx],ji,TJ-lagl)');
                jmx = jmx+1;
            endo;
            mmux = mmux.*.ones(TJ-lagl,1);
            mom = y-mu.*(my./mmu);
            db = -(mu.*(my./mmu)).*x
                   + (mu.*(my./(mmu^2))).*mmux;
        else;
            if addit==1;
                mom = y-mu;
                db = -mu.*x;
            else;
                mom = (y-mu)./mu;
                db = -(y./mu).*x;
            endif;
        endif;
    endif;
    retp(mom,z,db);
endp;


proc fas(bb,data);
    local a,b,dt,zmom,j,sji,ji,yx,
          mom,z,dab,f;


    if lfmy>0;
        a = bb[1:lfmlag];
        b = bb[lfmlag+1:rows(bb)];
    else;
        a = 0;
        b = bb;
    endif;
    open dt = ^dataset;
    zmom = 0;
    zgy = 0;
    j = 1;
    do until j>rows(datn);
        sji = 0;
        NJ = datn[j];
        TJ = datt[j];
        ji = minc(nind|NJ);
        do until sji==NJ;
            yx = readr(dt,ji*TJ);
            if TJ<lagl+llev+1;
                goto jpppp;
            endif;

            {mom,z,dab} = makemom(yx,ji,j,a,b);

            zmom = zmom + z'mom;
            zgy = zgy + z'dab;
jpppp:
            sji = sji+ji;
            ji = minc(nind|(NJ-sji));

        endo;
        j = j+1;
    endo;
    dt = close(dt);

    f = zmom'w*zmom;
    sco = zgy'w*zmom;
    hes = zgy'w*zgy;


    if f<staf;
        locate 8,1;
        "ITERATION";;
        mliter;
        print;
        j = 1;
        do until j>rows(bb);
            $varx[j];;bb[j];;
            if floor(j/3)==ceil(j/3);
                print;
            endif;
            j = j+1;
        endo;
        print;
        print;
        "FUNCTION VALUE";;
        format 12,4;print f;
        format 10,4;
        mliter = mliter+1;
        staf = f;
    endif;

    retp(-f);
endp;

proc fassco(b,data);
    retp(-2*sco');
endp;

proc fashes(b,data);
    retp(-2*hes);
endp;



proc (3)=like(abb,data);
    local iter,j,f,direc,step,ihes,btry,ftry;

    f = fas(abb,data);
    locate 1,1;

    __title;

    iter = 1;
    do while iter <=1000;
        locate 5,1;
        "ITERATION";;
        print iter;
        print;
        j = 1;
        do until j>rows(abb);
            $varx[j];;abb[j];;
            if floor(j/3)==ceil(j/3);
                print;
            endif;
            j = j+1;
        endo;
        print;
        print;
        "FUNCTION VALUE";;
        format 12,4;print -f;
        format 10,4;
        print;
        trap 1;
        ihes = invpd(hes);
        if scalerr(ihes);
            print "HESSIAN MATRIX IS NOT NEGATIVE DEFINITE";
            print "PROGRAMME ABORTED";
            end;
        endif;
        trap 0;
        direc = invpd(hes)*sco;
        if sco'direc <= 1e-10;
            retp(abb,f,iter);
        endif;
        step = 1;
        do while step >=.001;
            btry = abb - step*direc;
            ftry = fas(btry,data);
            if ftry <= f; goto nextstep; endif;
            abb = btry;
            f = ftry;
            goto stepfin;
nextstep:
            step = step/2;
        endo;
        print "THE ALGORITHM FAILED TO FIND A MINIMUM";
        print "PROGRAMME ABORTED";
        end;
        retp(0,0,0);
stepfin:
        iter = iter + 1;
    endo;

endp;

proc (2)=autotest(bb,avarb);
    local a,b,vv_1,vv_2,dabmom_1,dabmom_2,zmom1,zmom2,
          vv1,vv2,test1,test2,dt,j,sji,ji,yx,i,
          zi,momi,mom,z,dab,dab1,dab2,mom_1,mom_2,
          mom1,mom2,momi_1,momi_2,momi1,zi1,momi2,
          zi2,e,v1,v2,p;

    if lfmy>0;
        a = bb[1:lfmlag];
        b = bb[lfmlag+1:rows(bb)];
    else;
        a = 0;
        b = bb;
    endif;

    vv_1 = 0;
    vv_2 = 0;
    dabmom_1 = 0;
    dabmom_2 = 0;
    zmom1 = 0;
    zmom2 = 0;
    vv1 = 0;
    vv2 = 0;
    test1 = 0;
    test2 = 0;

    open dt = ^dataset;
    j = 1;
    do until j>rows(datn);
        sji = 0;
        NJ = datn[j];
        TJ = datt[j];
        ji = minc(nind|NJ);
        do until sji==NJ;
            yx = readr(dt,ji*TJ);
            if TJ<lagl+llev+1;
                goto jpppp;
            endif;

            {mom,z,dab} = makemom(yx,ji,j,a,b);

            if TJ>lagl+1;
                e = ones(rows(mom)/ji,1);
                e[1] = 0;
                e = ones(ji,1).*.e;
                mom_1 = selif(lag1(mom),e);
                mom1 = selif(mom,e);
                dab1 = selif(dab,e);
                vv_1 = vv_1 + mom_1'mom1;
                dabmom_1 = dabmom_1 + dab1'mom_1;
            endif;


            if TJ>lagl+2;
                e = ones(rows(mom)/ji,1);
                e[1] = 0;e[2]=0;
                e = ones(ji,1).*.e;
                mom_2 = selif(lagn(mom,2),e);
                mom2 = selif(mom,e);
                dab2 = selif(dab,e);
                vv_2 = vv_2 + mom_2'mom2;
                dabmom_2 = dabmom_2 + dab2'mom_2;
            endif;

            i = 1;
            do until i>ji;

                zi = z[(i-1)*(TJ-lagl)+1:i*(TJ-lagl),.];
                momi = mom[(i-1)*(TJ-lagl)+1:i*(TJ-lagl)];

                if TJ>lagl+1;
                    momi_1 = trimr(lag1(momi),1,0);
                    momi1 = trimr(momi,1,0);
                    zi1 = trimr(zi,1,0);
                    zmom1 = zmom1 + zi'momi*momi1'momi_1;
                    vv1 = vv1 + (momi_1'momi1)^2;
                endif;

                if TJ>lagl+2;
                    momi_2 = trimr(lagn(momi,2),2,0);
                    momi2 = trimr(momi,2,0);
                    zi2 = trimr(zi,2,0);
                    zmom2 = zmom2 + zi'momi*momi2'momi_2;
                    vv2 = vv2 + (momi_2'momi2)^2;
                endif;

                i = i+1;
            endo;


jpppp:
            sji = sji+ji;
            ji = minc(nind|(NJ-sji));

        endo;
        j = j+1;
    endo;
    dt = close(dt);


    if vv1>0;
        v1 = vv1 + 2*dabmom_1'invpd(hes)*zgy'w*zmom1 + dabmom_1'avarb*dabmom_1;
        test1 = vv_1/sqrt(v1);
    endif;
    if vv2>0;
        v2 = vv2 + 2*dabmom_2'invpd(hes)*zgy'w*zmom2 + dabmom_2'avarb*dabmom_2;
        test2 = vv_2/sqrt(v2);
    endif;

    p = 2*cdfnc(abs(test1|test2));

    retp(test1|test2,p);
endp;

if timevar==0;
    seqz = 0;
endif;


if nonseqz==0 and seqz==0 and ((timedum-timez)>-1);
    print "NO INSTRUMENTS SELECTED";
    print "PROGRAMME ABORTED";
    end;
endif;

if rows(lx)/=rows(xvar);
    print "# OF ELEMENTS IN lx NOT EQUAL TO # IN xvar";
    print "PROGRAMME ABORTED";
    end;
endif;

j = 1;
do until j>rows(xvar);
    if lx[j]<0;
        print "lx IS NEGATIVE for " $xvar[j];
        print "PROGRAMME ABORTED";
        end;
    endif;
    j = j+1;
endo;



if nonseqz>0;
    if rows(lnseqz)/=rows(nseqzvar);
        print "# OF ELEMENTS IN lnseqz NOT EQUAL TO # IN nseqzvar";
        print "PROGRAMME ABORTED";
        end;
    endif;
    j = 1;
    do until j>rows(nseqzvar);
        if lnseqz[j]<0;
            print "lx IS NEGATIVE for " $nseqzvar[j];
            print "PROGRAMME ABORTED";
            end;
        endif;
    j = j+1;
    endo;
endif;


if lfmy>0;
    if nonseqz>0;
        lagln = maxc(lx|lnseqz|lfmlag);
    else;
        lagln = maxc(lx|lfmlag);
    endif;
else;
    if nonseqz>0;
        lagln = maxc(lx|lnseqz);
    else;
        lagln = maxc(lx);
    endif;
endif;
if lagln>lagl;
    print;
    print "MODEL LAG LENGTH IS LARGER THAN lagl";
    print "lagl SET TO" lagln;
    lagl = lagln;
    pause(3);
endif;
if lagln<lagl;
    print "WARNING: lagl IS LARGER THAN MODEL LAG LENGTH";
    pause(3);
endif;


if model==2;
    lagl = lagl+1;
endif;


fake = ones(2,1);
let ss = "0" _1 _2 _3 _4 _5 _6 _7 _8 _9
         _10 _11 _12 _13 _14 _15 _16 _17 _18 _19
         _20 _21 _22 _23 _24 _25 _26 _27 _28 _29
         _30 _31 _32 _33 _34 _35 _36 _37 _38 _39
         _40 _41 _42 _43 _44 _45 _46 _47 _48 _49
         _50 _51 _52 _53 _54 _55 _56 _57 _58 _59
         _60 _61 _62 _63 _64 _65 _66 _67 _68 _69
         _70 _71 _72 _73 _74 _75 _76 _77 _78 _79
         _80 _81 _82 _83 _84 _85 _86 _87 _88 _89
         _90 _91 _92 _93 _94 _95 _96 _97 _98 _99;

sst = "1"|"2"|"3"|"4"|"5"|"6"|"7"|"8"|"9"|
         "10"|"11"|"12"|"13"|"14"|"15"|"16"|"17"|"18"|"19"|
         "20"|"21"|"22"|"23"|"24"|"25"|"26"|"27"|"28"|"29"|
         "30"|"31"|"32"|"33"|"34"|"35"|"36"|"37"|"38"|"39"|
         "40"|"41"|"42"|"43"|"44"|"45"|"46"|"47"|"48"|"49"|
         "50"|"51"|"52"|"53"|"54"|"55"|"56"|"57"|"58"|"59"|
         "60"|"61"|"62"|"63"|"64"|"65"|"66"|"67"|"68"|"69"|
         "70"|"71"|"72"|"73"|"74"|"75"|"76"|"77"|"78"|"79"|
         "80"|"81"|"82"|"83"|"84"|"85"|"86"|"87"|"88"|"89"|
         "90"|"91"|"92"|"93"|"94"|"95"|"96"|"97"|"98"|"99";



open dt = ^auxset;
datn = readr(dt,rowsf(dt));
dt = close(dt);

datt = datn[.,1];
datn = datn[.,2];
NTj = sumc(datt.*datn);
N = sumc(datn);
tmax = maxc(datt);

open dt = ^dataset;
nrdt = rowsf(dt);
dt = close(dt);

if NTj>nrdt;
    "THE NUMBER OF ROWS IN THE DATA SET IS SMALLER THAN";
    "INDICATED BY THE AUXILIARY FILE";
    "PROGRAMME ABORTED";
    END;
endif;

if NTj<nrdt;
    "THE NUMBER OF ROWS IN THE DATA SET IS LARGER THAN";
    "INDICATED BY THE AUXILIARY FILE";
    "PROGRAMME ABORTED";
    END;
endif;


j = 1;
do until j>rows(lx);
    if lx[j]>=tmax;
      print "lx IS AT LEAST AS LARGE AS THE LONGEST TIME SERIES FOR" $xvar[j];
        print "PROGRAMME ABORTED";
    endif;
    j = j+1;
endo;
if nonseqz>0;
    j = 1;
    do until j>rows(lnseqz);
        if lnseqz[j]>=tmax;
print "lnseqz IS AT LEAST AS LARGE AS THE LONGEST TIME SERIES FOR" $nseqzvar[j];
            print "PROGRAMME ABORTED";
        endif;
        j = j+1;
    endo;
endif;

if tmax<lagl+llev+1;
    print "THERE ARE NO OBERVATIONS THAT HAVE SUFFICIENT TIME SERIES";
    print "OBSERVATIONS FOR THE CURRENT lagl AND llev SETTINGS";
    print "PROGRAMME ABORTED";
    end;
endif;



{yvar,indy} = indices(dataset,yvar);
{xvar,indx} = indices(dataset,xvar);
uindx = indy|indx;

if seqz>0;
    if rows(lseqz1)/=rows(seqzvar);
        print "# OF ELEMENTS IN lseqz1 NOT EQUAL TO # IN seqzvar";
        print "PROGRAMME ABORTED";
        end;
    endif;
    if rows(lseqz2)/=rows(seqzvar);
        print "# OF ELEMENTS IN lseqz2 NOT EQUAL TO # IN seqzvar";
        print "PROGRAMME ABORTED";
        end;
    endif;
    {seqzvar,indz} = indices(dataset,seqzvar);
    uindx = uindx|indz;
    varz = seqzvar;
    j = 1;
    do until j>rows(varz);
        if lseqz1[j]>=tmax;
    print "lseqz1 IS AT LEAST AS LARGE AS THE LONGEST TIME SERIES FOR" $varz[j];
            print "PROGRAMME ABORTED";
            end;
        endif;
        if lseqz2[j]<lseqz1[j];
            print "lseqz2 IS SMALLER THAN lseqz1 FOR" $varz[j];
            print "PROGRAMME ABORTED";
            end;
        endif;
        if -lseqz1[j]>=tmax;
print "WARNING: -lseqz1 IS AT LEAST AS LARGE AS THE LONGEST TIME SERIES FOR"
$varz[j];
            print "-lseqz1 FOR THIS INSTRUMENTS IS SET TO" tmax-1;
            pause (3);
            lseqz1[j] = -(tmax-1);
        endif;

        if lseqz2[j]>=tmax;
print "WARNING: lseqz2 IS AT LEAST AS LARGE AS THE LONGEST TIME SERIES FOR"
$varz[j];
            print "lseqz2 FOR THIS INSTRUMENTS IS SET TO" tmax-1;
            pause (3);
            lseqz2[j] = tmax-1;
        endif;
        str = ""$+varz[j];
        sl = strlen(str);
        if sl>4;
            varz[j] = strsect(str,1,4);
        endif;
        if lseqz1[j]>=0;
            varz[j] = varz[j] $+ ss[lseqz1[j]+1];
        else;
            varz[j] = varz[j] $+ sst[-lseqz1[j]];
        endif;
        if lseqz2[j]>=0;
            varz[j] = varz[j] $+ ss[lseqz2[j]+1];
        else;
            varz[j] = varz[j] $+ sst[-lseqz2[j]];
        endif;
        j = j+1;
    endo;
endif;
if nonseqz>0;
    {nseqzvar,indiv} = indices(dataset,nseqzvar);
    uindx = uindx|indiv;
    variv = nseqzvar;
    j = 1;
    do until j>rows(variv);
        if lnseqz[j]>0;
            str = ""$+variv[j];
            sl = strlen(str);
            if sl>6;
                variv[j] = strsect(str,1,6);
            endif;
            variv[j] = variv[j] $+ ss[lnseqz[j]+1];
        elseif lnseqz[j]<0;
            str = ""$+variv[j];
            sl = strlen(str);
            if sl>6;
                variv[j] = strsect(str,1,6);
            endif;
            variv[j] = variv[j] $+ sst[-lnseqz[j]];
        endif;
        j = j+1;
    endo;
endif;

if seqz>0;
    if nonseqz>0;
        varz = varz|variv;
    endif;
elseif nonseqz>0;
    varz = variv;
endif;

if timevar/=0;
    {timevar,indyear} = indices(dataset,timevar);
    uindx = uindx|indyear;
endif;

uindx = uindx[uniqindx(uindx,1)];
__output = 0;
__miss = 2;
{uvar,mu,vu,stdu,minu,maxu,valu,misu} = dstat(dataset,uindx);
__output = 2;

if maxc(misu)>0;
    print "MISSING VALUES ENCOUNTERED FOR";
    print $selif(uvar,misu.>0)';
    print "PROGRAMME ABORTED";
    end;
endif;

dataseto = dataset;

if devvar==1;
    names = getname(dataset);
    uindx = indx[uniqindx(indx,1)];
    uxvar = names[uindx];
    xvaro = xvar;
    let dr = d;
    uxvar = dr $+ uxvar;
    xvar = dr $+ xvar;
    fname = "temp";
    names = names|uxvar;
    create dt = ^fname with ^names,0,8;
    dt = close(dt);

    __output = 0;
    {var,mx,vx,stdx,minx,maxx,valxt,misx} = dstat(dataset,0);
    __output = 2;


    open dtw = ^fname for append;
    open dtr = ^dataset;
    jd = 0;
    do until jd==rowsf(dtr);
       md = minc(100|(rowsf(dtr)-jd));
       data = readr(dtr,md);
       data = data~(data[.,uindx]-mx[uindx]');
       u = writer(dtw,data);
       jd = jd+md;
    endo;
    dtr = close(dtr);
    dtw = close(dtw);
    dataset = "temp";
    clearg data;
    {xvar,indx} = indices(dataset,xvar);
endif;


if sval>0;
    if rows(startvalx)==rows(xvar);
        b = startvalx;
    else;
        print;
        print "# OF STARTING VALUES NOT COMPATIBLE WITH # OF REGRESSORS";
        print "STARTING VALUES FOR REGRESSORS SET TO 0";
        pause(3);
        b = zeros(rows(xvar),1);
    endif;
else;
    b = zeros(rows(xvar),1);
endif;



if timevar/=0;
    __output=0;
   {timevar,myr,vyr,stdyr,minyear,maxyear,valyr,misyr} = dstat(dataset,timevar);
    __output=2;

    if ((maxyear-minyear-lagl+1)<1) or ((maxyear-minyear+1)<tmax);
        print $timevar " IS NOT PROPERLY SPECIFIED";
        print "PROGRAMME ABORTED";
        end;
    endif;
    years = makeyears;

else;
    years = 0;
    timedum = 0;
    timez = 0;
endif;

varx = xvar;
if devvar == 1;
    varx = xvaro;
endif;
j = 1;
do until j>rows(xvar);
    if lx[j]>0;
        str = ""$+varx[j];
        sl = strlen(str);
        len = 8-strlen(""$+ss[lx[j]+1]);
        if sl>len;
            varx[j] = strsect(str,1,len);
        endif;
        varx[j] = varx[j] $+ ss[lx[j]+1];
    endif;
    j = j+1;
endo;


if timedum==1;
    b = zeros(rows(years),1)|b;
    /* varx = ((ones(rows(years),1).*"T")$+sst[1+lagl:rows(years)+lagl])|varx;
       varz = ((ones(rows(years),1).*"T")$+sst[1+lagl:rows(years)+lagl])|varz; */
    varx = ((reshape(satocv("T"), rows(years), 1)) $+sst[1+lagl:rows(years)+lagl])|varx;
    varz = ((reshape(satocv("T"), rows(years), 1)) $+sst[1+lagl:rows(years)+lagl])|varz;
    if model==1;
        b = trimr(b,1,0);
        varx = trimr(varx,1,0);
        varz = trimr(varz,1,0);
    endif;
    if model==0;
        varx[1] = "CONST";
        b[1] = startvalc;
    endif;
elseif timez==1;
    /* varz = ((ones(rows(years),1).*"T")$+sst[1+lagl:rows(years)+lagl])|varz; */
    varz = ((reshape(satocv("T"), rows(years), 1)) $+sst[1+lagl:rows(years)+lagl])|varz;
    if model==1;
        varz = trimr(varz,1,0);
    endif;
endif;

if model==0 and timedum /= 1;
    if sval>0;
        b = startvalc|b;
    else;
        b = 0|b;
    endif;
    varx = "CONST"|varx;
    if timez /= 1;
        varz = "CONST"|varz;
    endif;
endif;

if lfmy>0;
    if sval>0;
        if rows(startvaly)==lfmlag;
            b = startvaly|b;
        else;
            print;
            "# OF STARTING VALUES FOR ";;
            "LAGGED DEP. VAR. NOT COMPATIBLE WITH lfmlag";
            print "STARTING VALUES FOR LAGGED DEP. VAR. SET TO 0";
            pause(3);
            b = zeros(lfmlag,1)|b;
        endif;
    else;
        b = zeros(lfmlag,1)|b;
    endif;

    str = ""$+yvar;
    sl = strlen(str);
    len = 8-strlen(""$+ss[lfmlag+1]);
    vary = strsect(str,1,len);

    /* vary = ones(lfmlag,1)*vary; */
    vary = reshape(satocv(vary), lfmlag, 1);
    vary = vary $+ ss[2:lfmlag+1];
    varx = vary[1:lfmlag]|varx;
endif;

_max_Parnames = varx;



trap 1;
w = invpd(makezz);
if scalerr(w);
    print;
    print "THE WEIGHT MATRIX Z'Z IS SINGULAR";
    print "PROGRAMME ABORTED";
    end;
endif;
trap 0;

if cols(w)<rows(b);
    print;
    print "# OF INSTRUMENTS SMALLER THAN # OF PARAMETERS";
    print "PROGRAMME ABORTED";
    end;
elseif cols(w)==rows(b);
    justid = 1;
endif;


cls;

__title = "ONE-STEP ESTIMATION";
locate 1,1;
__title;
locate 3,1;
if model==0;
    "LEVELS MODEL ";;
    if addit==1;
        "ADDITIVE MOMENT CONDITIONS ";
    else;
        "MULTIPLICATIVE MOMENT CONDITIONS";
    endif;
elseif model==1;
    "WITHIN GROUP MEAN SCALING";
else;
    "QUASI DIFFERENCING ";;
    if qdif==0;
        "WOOLDRIDGE MOMENT CONDITIONS";
    else;
        "CHAMBERLAIN MOMENT CONDITIONS";
    endif;
    if devvar>0;
        "EXPLANATORY VARIABLES IN DEVIATION OF OVERALL MEANS";
    endif;

endif;
print;

{b1,f1,iter1} = like(b,fake);

trap 1;
ihes = invpd(hes);
if scalerr(ihes);
    print;
    print "THE HESSIAN MATRIX OF THE ONE-STEP ESTIMATOR IS SINGULAR";
    if model==2 and qdif==1;
        print "SETTING DEVVAR=1 MAY HELP";
    endif;
    print "PROGRAMME ABORTED";
    end;
endif;
trap 0;

w2 = makezuz(b1);
avarb1 = ihes*(zgy'w*w2*w*zgy)*ihes;
seb1 = sqrt(diag(avarb1));
pb1 = 2*cdfnc(abs(b1./seb1));
if (model/=1) and (timevar/=0);
    {act1,p1} = autotest(b1,avarb1);
endif;


trap 1;
w = invpd(w2);
if scalerr(w);
    print;
    print "THE EFFICIENT WEIGHT MATRIX Z'uu'Z IS SINGULAR";
    print "SETTING DEVVAR=1 WITH WOOLDRIDGE MOMENT CONDITIONS MAY HELP";
    print "PROGRAMME ABORTED";
    end;
endif;
trap 0;


if justid /= 1;

    __title = "TWO-STEP ESTIMATION";
    locate 1,1;
    __title;
    {b2,f2,iter2} = like(b1,fake);


    trap 1;
    ihes = invpd(hes);
    if scalerr(ihes);
        print;
        print "THE HESSIAN MATRIX OF THE TWO-STEP ESTIMATOR IS SINGULAR";
        print "SETTING DEVVAR=1 WITH WOOLDRIDGE MOMENT CONDITIONS MAY HELP";
        print "PROGRAMME ABORTED";
    end;
    endif;
    trap 0;


    avarb2 = ihes;
    seb2 = sqrt(diag(avarb2));
    pb2 = 2*cdfnc(abs(b2./seb2));
    if (model/=1) and (timevar/=0);
        {act2,p2} = autotest(b2,avarb2);
    endif;
endif;

let varm = m1 m2;


let mask[1,3] = 0 1 1;
let fmt[3,3] = "-*.*s" 12 8  "*.*lf" 10 4 "*.*lf" 10 4;

laglo = lagl;
if model==2;
    laglo = laglo-1;
endif;

if devvar==1;
    xvar = xvaro;
endif;

dataset = dataseto;

pdate = date;

format /rd 8,0;

output on;

print;
"/**** EXPEND                                                           ****/ ";
"/**** nonlinear gmm for EXPonential models with ENDogenous regressors  ****/ ";
"/**** Version 1.02, July 2003                                          ****/ ";
"/**** Frank Windmeijer, Institute for Fiscal Studies                   ****/ ";
"/**** f.windmeijer@ifs.org.uk                                          ****/ ";
print;
"----------------------------------------------------------------------------";
if model==1;
"WITHIN GROUP MEAN SCALING ESTIMATION";
elseif model==2;
    if qdif==1;
"CHAMBERLAIN MOMENT CONDITIONS";
    else;
"WOOLDRIDGE MOMENT CONDITIONS";
    endif;
else;
"MODEL IN LEVELS";
    if addit==1;
"ADDITIVE MOMENT CONDITIONS";
    else;
"MULTIPLICATIVE MOMENT CONDITIONS";
    endif;
endif;
"DATASET        " $dataset;
if devvar>0;
"EXPLANATORY VARIABLES IN DEVIATION OF OVERALL MEANS";
endif;
"DATE           " datestr(0);
"N              " sumc(selif(datn,datt.>lagl+llev)) "  NT             "
sumc(selif(datn,datt.>lagl+llev).*(selif(datt,datt.>lagl+llev)-lagl));
"LAGL           " laglo "  LLEV           " llev;
if timevar/=0;
"PERIOD         " (minyear+lagl)~maxyear;
endif;
"DEP. VAR.      " $yvar;
"INSTRUMENTS    ";;
format /ldn;
if timedum==1 or timez==1; format 4,0; $varz[1:rows(years)-(model==1)]';
"               ";;
elseif model==0; format 6,0; $varz[1];"               ";;endif;
if seqz>0;
    j = 1;
    do until j>rows(seqzvar);
        str = ""$+seqzvar[j];
        sl = strlen(str);
        format sl,0;
        $seqzvar[j];;
        if lseqz1[j]>=0;
            format strlen(""$+ss[lseqz1[j]+1]),0;
            $ss[lseqz1[j]+1];;
        else;
            format strlen(""$+sst[-lseqz1[j]+1]),0;
            $sst[-lseqz1[j]];;
        endif;
        if lseqz2[j]>=0;
            format strlen(""$+ss[lseqz2[j]+1])+1,0;
            $ss[lseqz2[j]+1];;
        else;
            format strlen(""$+sst[-lseqz2[j]+1])+1,0;
            $sst[-lseqz2[j]];;
        endif;
        j = j+1;
        if floor(j/5)==ceil(j/5);print;"               ";;endif;
    endo;
    if (floor(j/5))/=(ceil(j/5));print;"               ";;endif;
endif;
if nonseqz>0;
    j = 1;
    do until j>rows(nseqzvar);
        str = ""$+nseqzvar[j];
        sl = strlen(str);
        if lnseqz[j]==0;
            format sl+1,0;
        else;
            format sl,0;
        endif;
        $nseqzvar[j];;
        if lnseqz[j]>0;
            format strlen(""$+ss[lnseqz[j]+1])+1,0;
            $ss[lnseqz[j]+1];;
        elseif lnseqz[j]<0;
            format strlen(""$+sst[-lnseqz[j]+1])+1,0;
            $sst[-lnseqz[j]];;
        endif;
        j = j+1;
        if floor(j/5)==ceil(j/5);print;"               ";;endif;
    endo;
endif;
print;

format /rd;

format 8,4;
if justid /= 1;
    "SARGAN DOF P:  " (-f2);;
    format 4,0; (cols(w)-rows(b2));;
    format 8,4; cdfchic(-f2,cols(w)-rows(b2));
endif;
"----------------------------------------------------------------------------";
"                      ONE-STEP ";
"# ITERATIONS";;
if iter1<10;
    format 2,0;
elseif iter1<100;
    format 3,0;
elseif iter1<1000;
    format 4,0;
elseif iter1<10000;
    format 5,0;
else;
    format 8,0;
endif;
iter1;print;

"                 coeff    rob se   t-ratio   p-value ";

if lfmy>0;
   ylen = strlen(""$+yvar);
   j = 1;
   do until j>lfmlag;
      format /ldn ylen,0;
      $yvar;;
      format 12-ylen,0;$ss[j+1];;
      format /rdn 10,4;
      b1[j];;seb1[j];;(b1[j]/seb1[j]);;pb1[j];
      j = j+1;
   endo;
endif;
if lfmy==0;
   j = 1;
endif;
if timedum==1;
   jj = 1;
   do until jj>rows(years)-(model==1);
      format /ldn 12,0;
      $varx[j];;
      format /rdn 10,4;
      b1[j];;seb1[j];;(b1[j]/seb1[j]);;pb1[j];
      j = j+1;
      jj = jj+1;
   endo;
elseif model==0;
   format /ldn 12,0;
   $varx[j];;
   format /rdn 10,4;
   b1[j];;seb1[j];;(b1[j]/seb1[j]);;pb1[j];
   j = j+1;
endif;
jj = 1;
do until j>rows(varx);
   format /ldn;
   if lx[jj]==0;
      format 12,0;
   else;
      xlen = strlen(""$+xvar[jj]);
      format xlen,0;
   endif;
   $xvar[jj];;
   if lx[jj]>0;
      format 12-xlen,0;
      $ss[lx[jj]+1];;
   endif;
   format /rdn;
   format 10,4;
   b1[j];;seb1[j];;(b1[j]/seb1[j]);;pb1[j];
   j = j+1;
   jj = jj+1;
endo;

print;


if (model/=1) and (timevar/=0);
"   TESTS FOR SERIAL CORRELATION AND P-VALUES ";
d=printfm(varm~act1~p1,mask,fmt);print;
endif;
print;
if justid == 1;
"MODEL JUST IDENTIFIED, ONE-STEP AND TWO-STEP COINCIDE" ;
else;
"                      TWO-STEP ";
"# ITERATIONS";;
if iter2<10;
    format 2,0;
elseif iter2<100;
    format 3,0;
elseif iter2<1000;
    format 4,0;
elseif iter2<10000;
    format 5,0;
else;
    format 8,0;
endif;
iter2;print;

"                 coeff    rob se   t-ratio   p-value ";




if lfmy>0;
   ylen = strlen(""$+yvar);
   j = 1;
   do until j>lfmlag;
      format /ldn ylen,0;
      $yvar;;
      format 12-ylen,0;$ss[j+1];;
      format /rdn 10,4;
      b2[j];;seb2[j];;(b2[j]/seb2[j]);;pb2[j];
      j = j+1;
   endo;
endif;
if lfmy==0;
   j = 1;
endif;
if timedum==1;
   jj = 1;
   do until jj>rows(years)-(model==1);
      format /ldn 12,0;
      $varx[j];;
      format /rdn 10,4;
      b2[j];;seb2[j];;(b2[j]/seb2[j]);;pb2[j];
      j = j+1;
      jj = jj+1;
   endo;
elseif model==0;
   format /ldn 12,0;
   $varx[j];;
   format /rdn 10,4;
   b2[j];;seb2[j];;(b2[j]/seb2[j]);;pb2[j];
   j = j+1;
endif;
jj = 1;
do until j>rows(varx);
   format /ldn;
   if lx[jj]==0;
      format 12,0;
   else;
      xlen = strlen(""$+xvar[jj]);
      format xlen,0;
   endif;
   $xvar[jj];;
   if lx[jj]>0;
      format 12-xlen,0;
      $ss[lx[jj]+1];;
   endif;
   format /rdn;
   format 10,4;
   b2[j];;seb2[j];;(b2[j]/seb2[j]);;pb2[j];
   j = j+1;
   jj = jj+1;
endo;

print;
if (model/=1) and (timevar/=0);
"   TESTS FOR SERIAL CORRELATION AND P-VALUES ";
d=printfm(varm~act2~p2,mask,fmt);print;
endif;
endif;

format 8,2;
print " Execution time is " (hsec-speed)/6000 " minutes" ;

if saveres==1;
    save b1;
    save v1 = avarb1;
    if justid /= 1;
        save b2;
        save v2 = avarb2;
    endif;
endif;


output off;



