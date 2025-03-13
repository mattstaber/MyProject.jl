module NPCG

include("AuxDynamics.jl")
include("BankReversal.jl")
include("EntryAnalysis.jl")
include("calAirdens.jl") 
include("calDelpsi.jl")
include("calDrdc.jl")
include("calJ.jl") 
include("calPci2pcr.jl") 
include("calSphdis.jl") 
include("EntryAnalysis.jl")


function npcg(tc, ec, TU, GAT, DU, siglmt, tstp, sigdlmt, sigddlmt, sigcmdprv, sigcmdprv2, BRprev, sigs, BAL, BRL, X0, thetatgt, phitgt, Case_Number, auxdata)
    if tc <= GAT / TU
        return 0, sigs, BRprev, sigcmdprv, sigcmdprv2
    else
        inum = 0
        
        while true
            G = cal_z(X0, sigs, ec, auxdata)
            
            if abs(G) <= 1e3 / DU
                break
            end
            
            J = cal_J(X0, sigs, ec, auxdata)
            
            sigs -= G / J
            
            inum += 1
            println("   Correction at t = ", round(tc * TU), " s")
            
            if norm(J) < 1e-6
                break
            end
        end
                
        if BAL == 1
            if abs(sigs) >= siglmt
                sigs = siglmt
            end
        end
        
        if BRL == 1
            BR = cal_BR(X0, BRprev, Case_Number, auxdata)
        elseif BRL == 2
            BR = cal_BR_prdt(X0, ec, sigs, BRprev, auxdata)
        end
        
        if BR != BRprev
            Rgo = DU * 1e-3 * cal_sphdis(X0[2], X0[3], thetatgt, phitgt)
            println("       Bank reversal at Rgo = ", round(Rgo), " km")
        end
        
        BRprev = BR
        bankcmd = sigs[1] * BR
        
        if BAL == 2
            minbnk = max(sigcmdprv - sigdlmt * tstp, 2 * sigcmdprv - sigcmdprv2 - sigddlmt * tstp^2)
            maxbnk = min(sigcmdprv + sigdlmt * tstp, 2 * sigcmdprv - sigcmdprv2 + sigddlmt * tstp^2)
            
            bankcmd = clamp(bankcmd, minbnk, maxbnk)
            
            if abs(bankcmd) >= siglmt
                bankcmd = siglmt
            end
            
            sigcmdprv2 = sigcmdprv
            sigcmdprv = bankcmd
        end
        
        return bankcmd, sigs, BRprev, sigcmdprv, sigcmdprv2
    end
end

end