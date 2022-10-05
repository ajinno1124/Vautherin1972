module NuclParameters
    mutable struct NuclParams
        t0::Float64
        t1::Float64
        t2::Float64
        t3::Float64
        x0::Float64
        x1::Float64
        x2::Float64
        x3::Float64
        σ::Float64
        W0::Float64
    end

    function getParams(ParamType::String)
        #SLy4
        if ParamType=="SLy4"
            t0=-2488.91
            t1=486.82
            t2=-546.39
            t3=13777.0
            x0=0.834
            x1=-0.344
            x2=-1.000
            x3=1.354
            σ=1/6
            W0=123.0
        elseif ParamType=="SLy5"
            t0=-2484.88
            t1=483.13
            t2=-549.40
            t3=13763.0
            x0=0.778
            x1=-0.328
            x2=-1.000
            x3=1.267
            σ=1/6
            W0=126.0
        elseif ParamType=="SLy6"
            t0=-2479.50
            t1=462.18
            t2=-448.61
            t3=13763.0
            x0=0.825
            x1=-0.465
            x2=-1.000
            x3=1.355
            σ=1/6
            W0=122.0
        elseif ParamType=="SLy7"
            t0=-2482.41
            t1=457.97
            t2=-419.85
            t3=13677.0
            x0=0.846
            x1=-0.511
            x2=-1.000
            x3=1.391
            σ=1/6
            W0=126.0
        elseif ParamType=="SkM*"
            t0=-2645.00
            t1=410.00
            t2=-135.00
            t3=15595.0
            x0=0.09
            x1=0.00
            x2=0.00
            x3=0.00
            σ=1/6
            W0=130.0
        elseif ParamType=="VB1"
            t0=-1057.3
            t1=235.9
            t2=-100.0
            t3=14463.5
            x0=0.56
            x1=0.00
            x2=0.00
            x3=0.00
            σ=1
            W0=120.0
        elseif ParamType=="SKS3"
            t0=-2014.7
            t1=361.0
            t2=-29.5
            t3=12756
            x0=-0.319
            x1=0.732
            x2=4.95
            x3=-0.904
            σ=0.2604
            W0=94
        end

        return NuclParams(t0,t1,t2,t3,x0,x1,x2,x3,σ,W0)
    end

    function getaN(ParamType::String)
        aN=zeros(Float64,10)
        p=getParams(ParamType)
        aN[1]=0.25*p.t0*(2+p.x0)
        aN[2]=-0.25*p.t0*(2*p.x0+1)
        aN[3]=1/24*p.t3*(2+p.x3)
        aN[4]=-1/24*p.t3*(2*p.x3+1)
        aN[5]=1/8*(p.t1*(2+p.x1)+p.t2*(2+p.x2))
        aN[6]=1/8*(p.t2*(2*p.x2+1)-p.t1*(2*p.x1+1))
        aN[7]=1/32*(3*p.t1*(2+p.x1)-p.t2*(2+p.x2))
        aN[8]=-1/32*(3*p.t1*(2*p.x1+1)+p.t2*(2*p.x2+1))
        aN[9]=-1/16*(p.t1*p.x1+p.t2*p.x2)
        aN[10]=1/16*(p.t1-p.t2)

        if ParamType=="VB1"
            aN[3]=p.t3/8
            aN[4]=-aN[3]
        end

        return aN
    end
end


module LambdaParameters
    mutable struct LambdaParams
        γ::Float64
        u0::Float64
        u1::Float64
        u2::Float64
        u3::Float64
        u3p::Float64
        y0::Float64
        y3::Float64
    end

    function getParams(ParamType::String)
        if ParamType=="HPL1"
            γ = 1
            u0 = -326.395
            u1 = 72.627
            u2 = -8.584
            u3 = 0.0
            u3p = 1746.041
            y0 = -0.223
            y3 = -0.389
        elseif ParamType=="HPL2"
            γ = 1
            u0 = -399.946
            u1 = 83.426
            u2 = 11.455
            u3 = 0.0
            u3p = 2046.818
            y0 = -0.486
            y3 = -0.660
        elseif ParamType=="HPL3"
            γ =  1/3
            u0 = -498.515
            u1 = 65.203
            u2 = 19.001
            u3 = 0.0
            u3p = 995.832
            y0 = -0.492
            y3 = -0.444
        elseif ParamType=="HPL4"
            γ = 1/3
            u0 = -475.584
            u1 = 99.058
            u2 = -20.890
            u3 = 0.0
            u3p = 1375.172
            y0 = -0.350
            y3 = -0.724
        elseif ParamType=="NL1"
            γ = 1
            u0 = -253.3250
            u1 = 147.1264
            u2 = -83.5843
            u3 = 0.0
            u3p = 1684.9876
            y0 = 0.5802
            y3 = 0.4831
        elseif ParamType=="OL1"
            γ = 1
            u0 = -236.5835
            u1 = 116.8704
            u2 = -112.8812
            u3 = 0.0
            u3p = 1453.3493
            y0 = 0.1271
            y3 = -0.3110
        elseif ParamType=="NL2"
            γ = 1/3
            u0 = -518.620
            u1 = 82.0944
            u2 = -19.9772
            u3 = 0.0
            u3p = 1190.1894
            y0 = -0.1392
            y3 = 0.3126
        elseif ParamType=="OL2"
            γ = 1/3
            u0 = -417.7593
            u1 = 1.5460
            u2 = -3.2617
            u3 = 0.0
            u3p = 1102.2221
            y0 = -0.3854
            y3 = -0.5645
        elseif ParamType=="SKSH1"
            γ = 0
            u0 = -176.5
            u1 = -35.8
            u2 = 44.1
            u3 = 0.0
            u3p = 0.0
            y0 = 0.0
            y3 = 0.0
        elseif ParamType=="SKSH2"
            γ = 0
            u0 = -290.0
            u1 = 21.7
            u2 = -20.3
            u3 = 1850
            u3p = 0.0
            y0 = 0.0
            y3 = 0.0
        elseif ParamType=="GKW2"
            γ = 1/3
            u0 = -968.125
            u1 = -125.44938856118907
            u2 = -376.34816568356723
            u3 = 0.0
            u3p = 4371.717378386592
            y0 = 0.0
            y3 = 0.0
        elseif ParamType=="GKW3"
            γ = 1/3
            u0 = -500.62499999999994
            u1 = 295.4508964244826
            u2 = 886.3526892734478
            u3 = 0.0
            u3p = 4.912041998187181
            y0 = 0.0
            y3 = 0.0
        end

        return LambdaParams(γ,u0,u1,u2,u3,u3p,y0,y3)
    end

    function getaL(ParamType::String)
        aL=zeros(Float64,5)
        p=getParams(ParamType)
        aL[1]=p.u0*(1+0.5*p.y0)
        aL[2]=0.25*(p.u1+p.u2)
        aL[3]=1.0/8.0*(3*p.u1-p.u2)
        #if ParamType=="SKSH1" || ParamType=="SKSH2"
        #    aL[3]=0.25*(3*p.u1-p.u2)
        #end
        aL[4]=3.0/8.0*p.u3p*(1+0.5*p.y3)
        aL[5]=0.25*p.u3

        return aL
    end

end
