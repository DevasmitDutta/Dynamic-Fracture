### A Pluto.jl notebook ###
# v0.19.22

using Markdown
using InteractiveUtils

# ╔═╡ 9892a47e-b7fb-11ed-1236-2f700b706c6b
begin 
   using Pkg
   Pkg.activate("../")
   Pkg.instantiate()
   Pkg.status()
end

# ╔═╡ 601d363c-9fdd-4093-9fa5-b88ba7ce6cba
begin 
   using Gmsh: gmsh
   using GridapGmsh
   using Gridap
   using Gridap.Geometry
   using Gridap.TensorValues
   using Plots
   using LinearAlgebra
   using Gridap.ReferenceFEs
   using Gridap.Arrays
end

# ╔═╡ 35750bbc-d306-4806-9de0-0a3bbdaebcdb
begin
   using PlutoUI
   PlutoUI.TableOfContents()

end

# ╔═╡ a2764d82-0497-4111-818e-9a0ffbb62796
md"""
## Activating the Project

"""

# ╔═╡ 7a3de35c-6ba9-48ce-bbfb-e761a60cbbde
md"""
## Creating the mesh

"""

# ╔═╡ e36ffeea-13fa-431f-ae48-b44a6fedc76c
md"""
## Geometrical constraints of our model
$(Resource("https://imgur.com/UXruMa7.png", :width=> 600))
"""

# ╔═╡ b897e239-36ed-443a-beb8-965d46b1c0e5
const ls = 0.8

# ╔═╡ 9a44e0aa-e92d-486d-920e-93918ae0f791
begin
	const L = 228
	const H = 76.2
	const γₜ = 0.77
	const γ = 0.6 # [0.82, 0.6]
	const L1 = 101.6
	const L2 = γ*L1
	const CP = 19.4 # crack height
	const CL = (L/2)-L2 # crack location along length direction
	const CH = L/1000
	const hfc = ls/50
	const hf = ls/2.1
	const h = 3*hf
	const FMR = 5*ls
	
	gmsh.initialize()
	gmsh.option.setNumber("General.Terminal", 1)
	p1 = gmsh.model.geo.addPoint(0.0, 0.0, 0.0, h)
	p2 = gmsh.model.geo.addPoint((L/2)-L1-h, 0.0, 0.0, h)
	p3 = gmsh.model.geo.addPoint((L/2)-L1+h, 0.0, 0.0, h)
	p4 = gmsh.model.geo.addPoint(CL-CH, 0.0, 0.0, hf) 
	p5 = gmsh.model.geo.addPoint(CL-CH, CP, 0.0, hfc) 
	p6 = gmsh.model.geo.addPoint(CL+CH, CP, 0.0, hfc) 
	p7 = gmsh.model.geo.addPoint(CL+CH, 0.0, 0.0, hf)
	p8 = gmsh.model.geo.addPoint((L/2)+L1-h, 0.0, 0.0, h)
	p9 = gmsh.model.geo.addPoint((L/2)+L1+h,0.0, 0.0, h)
	p10 = gmsh.model.geo.addPoint(L, 0, 0.0, h)
	p11 = gmsh.model.geo.addPoint(L, H, 0.0, h)
	p12 = gmsh.model.geo.addPoint((L/2)+2*h, H, 0.0, hf)
	p13 = gmsh.model.geo.addPoint((L/2)-2*h, H, 0.0, hf)
	p14 = gmsh.model.geo.addPoint(0, H, 0.0, h)
	
	l1 = gmsh.model.geo.addLine(p1, p2)
	l2 = gmsh.model.geo.addLine(p2, p3)
	l3 = gmsh.model.geo.addLine(p3, p4)
	l4 = gmsh.model.geo.addLine(p4, p5)
	l5 = gmsh.model.geo.addLine(p5, p6)
	l6 = gmsh.model.geo.addLine(p6, p7)
	l7 = gmsh.model.geo.addLine(p7, p8)
	l8 = gmsh.model.geo.addLine(p8, p9)
	l9 = gmsh.model.geo.addLine(p9, p10)
	l10 = gmsh.model.geo.addLine(p10, p11)
	l11 = gmsh.model.geo.addLine(p11, p12)
	l12 = gmsh.model.geo.addLine(p12, p13)
	l13 = gmsh.model.geo.addLine(p13, p14)
	l14 = gmsh.model.geo.addLine(p14, p1)
	
	cl = gmsh.model.geo.addCurveLoop([l1,l2,l3,l4,l5,l6,l7,l8,l9,l10,l11,l12,l13,l14])
	ps = gmsh.model.geo.addPlaneSurface([cl])
	pg1 = gmsh.model.addPhysicalGroup(2, [ps])
	pg2 = gmsh.model.addPhysicalGroup(1, [l2]) 
	pg3 = gmsh.model.addPhysicalGroup(1,[l8])
	pg4 = gmsh.model.addPhysicalGroup(1,[l12])
	
	gmsh.model.setPhysicalName(2, pg1, "Domain")
	gmsh.model.setPhysicalName(1, pg2, "LeftSupp")
	gmsh.model.setPhysicalName(1, pg3, "RightSupp")
	gmsh.model.setPhysicalName(1, pg4, "LoadLine")
	
	p15 = gmsh.model.geo.addPoint(L/2, H, 0.0, hf)
	p16 = gmsh.model.geo.addPoint(L/2, 0.0, 0.0, hf)
	l16 = gmsh.model.geo.addLine(p6, p15)
	l17 = gmsh.model.geo.addLine(p15, p16)
	
	gmsh.model.mesh.field.add("Distance", 1)
	gmsh.model.mesh.field.setNumbers(1, "EdgesList", [l16,l17])
	
	gmsh.model.mesh.field.add("Threshold", 2)
	gmsh.model.mesh.field.setNumber(2, "IField", 1)
	gmsh.model.mesh.field.setNumber(2, "LcMin", hf)
	gmsh.model.mesh.field.setNumber(2, "LcMax", h)
	gmsh.model.mesh.field.setNumber(2, "DistMin", 0.15*H)
	gmsh.model.mesh.field.setNumber(2, "DistMax", 0.2*H)
	
	gmsh.model.mesh.field.setAsBackgroundMesh(2)
	
	gmsh.model.geo.synchronize()
	gmsh.model.mesh.generate(2)
	gmsh.write("JohnShah.msh")
	gmsh.finalize()
end

# ╔═╡ 08d10d55-bc28-4ac6-afb7-1d079a62548d
begin
	model = GmshDiscreteModel("JohnShah.msh")
    writevtk(model ,"JohnShah")
end

# ╔═╡ 13a77d13-7285-497f-a466-e0ac1685d9a5
md"""
## Defining the material parameters

"""

# ╔═╡ d39cc524-02b7-430d-ac2a-b0e2ffbe9e41
begin 
    const E_mat = 31.37*1e3 
    const ν_mat = 0.2   
    const ρ = 2400*1e-9

    G₁₂_mat = E_mat/(2*(1+ν_mat))
    const μ = G₁₂_mat
    const λ = (E_mat*ν_mat)/((1+ν_mat)*(1-2*ν_mat))
    const k = λ + μ

    J = 0.0
end

# ╔═╡ af71ff02-1d8b-470a-bbf2-72174eeadeaa
begin
	const Gc = 31.1*1e-3 
    const η = 1e-15
end

# ╔═╡ 9621fd3c-15a3-4403-b805-f49245a16726
begin
	function ElasFourthOrderConstTensor(E₁,E₂,ν₁₂,ν₂₁,G₁₂,PlanarState)
	    # 1 for Plane Stress and 2 Plane Strain Condition 
	  if PlanarState == 1
	      C1111 = E₁/(1-ν₁₂*ν₂₁)
	      C1122 = (ν₂₁*E₁)/(1-ν₁₂*ν₂₁)
	      C1112 = 0.0
	      C2222 = E₂/(1-ν₁₂*ν₂₁)
	      C2212 = 0.0
	      C1212 = G₁₂     
	  elseif PlanarState == 2
	      C1111 = (E₁*(1-ν₁₂*ν₂₁))/((1+ν₁₂)*(1-ν₁₂-2*ν₁₂*ν₂₁))
	      C1122 = (ν₂₁*E₁)/(1-ν₁₂-2*ν₁₂*ν₂₁)
	      C1112 = 0.0
	      C2222 = (E₂*(1-ν₁₂))/(1-ν₁₂-2*ν₁₂*ν₂₁)
	      C2212 = 0.0
	      C1212 = G₁₂ 
	  end
	  C_ten = SymFourthOrderTensorValue(C1111,C1112,C1122,C1112,C1212,C2212,C1122,C2212,C2222)
	    return  C_ten
	end
	const C_mat = ElasFourthOrderConstTensor(E_mat,E_mat,ν_mat,ν_mat,G₁₂_mat,2)
end

# ╔═╡ 1a22a2ba-3074-4b3a-9229-c0195645177a
md"""
## Defining the spectral-decomposition method

"""

# ╔═╡ 0af2040c-b6c4-42f3-848c-2def4b3e6881
### Spectral

function σ_mod(ε, ε_in, s_in)
            σM = (s_in^2 + η)*(λ*tr(ε)*one(ε) + 2*μ*(ε))
    return σM
end

# ╔═╡ 3d5dad6b-1db1-4dd6-b01f-23a5de21f9b0
### Spectral

function Eigen(ε)
    εArray = get_array(ε)
    Λ, P = eigen(εArray)
    ε1 = Λ[1]
    ε2 = Λ[2]
    if ε1 >= 0 &&  ε2 >= 0
        Λpos = [ε1 0; 0 ε2]
        Λneg = [0 0; 0 0]
    elseif ε1 >= 0 &&  ε2 < 0
        Λpos = [ε1 0; 0 0]
        Λneg = [0 0; 0 ε2]
    elseif ε1 < 0 &&  ε2 >= 0
        Λpos = [0 0; 0 ε2]
        Λneg = [ε1 0; 0 0]
    elseif ε1 < 0 &&  ε2 < 0
        Λpos = [0 0; 0 0]
        Λneg = [ε1 0; 0 ε2]
    end 
   εPos = P*Λpos*P'
   εNeg = P*Λneg*P'
   EpPos = TensorValue(εPos)
   EpNeg = TensorValue(εNeg)    
    return  EpPos, EpNeg
end


# ╔═╡ 753a4974-059e-432f-b344-13504b9d401c
function ψPos(ε_in)
    ϵ_Pos, ϵ_Neg = Eigen(ε_in)
    if tr(ε_in) >= 0 
        ψPlus = 0.5*λ*(tr(ε_in))^2 + μ*(ϵ_Pos ⊙ ϵ_Pos) 
    elseif tr(ε_in) < 0
        ψPlus = μ*(ϵ_Pos ⊙ ϵ_Pos)
    end
    return ψPlus
end

# ╔═╡ 1467f800-5ea9-4724-8749-81f2c4c467e7
function new_EnergyState(ψPlusPrev_in,ψhPos_in)
  ψPlus_in = ψhPos_in
  if ψPlus_in >= ψPlusPrev_in
    ψPlus_out = ψPlus_in
  else
    ψPlus_out = ψPlusPrev_in
  end
  true,ψPlus_out
end

# ╔═╡ 2a0d35c5-66e2-408c-b87d-f195eff673f6
function project(q,model,dΩ,order)
  reffe = ReferenceFE(lagrangian,Float64,order)
  V = FESpace(model,reffe,conformity=:L2)
  a(u,v) = ∫( u*v )*dΩ
  l(v) = ∫( v*q )*dΩ
  op = AffineFEOperator(a,l,V,V)
  qh = solve(op)
  qh
end

# ╔═╡ d8265d0e-a3e9-4a7f-8792-9275222b1810
begin
	order = 1
	degree = 2*order
end

# ╔═╡ 790de406-a97e-4988-896d-e8d62115cf36
begin
	Ω = Triangulation(model)
	dΩ = Measure(Ω,degree)
end

# ╔═╡ 830e4d8e-8f22-4eed-9e0d-d0c3a943a72c
begin
        labels = get_face_labeling(model)
        LoadTagId = get_tag_from_name(labels,"LoadLine")
        Γ_Load = BoundaryTriangulation(model,tags = LoadTagId)
        dΓ_Load = Measure(Γ_Load,degree)
        n_Γ_Load = get_normal_vector(Γ_Load)
end

# ╔═╡ b1d66ae6-5414-424e-b18d-c503082df7a8
md"""
## Defining the FE-Space for phase-field variable

"""

# ╔═╡ e8707054-6662-48a1-b758-3438990a5a48
begin
	reffe_PF = ReferenceFE(lagrangian,Float64,order)
	V0_PF = TestFESpace(model,reffe_PF;
	  conformity=:H1)
	U_PF = TrialFESpace(V0_PF)
end

# ╔═╡ d2b2cb89-2d65-4a27-b34d-b8576218198e
begin
	reffe_Disp = ReferenceFE(lagrangian,VectorValue{2,Float64},order)
        V0_Disp = TestFESpace(model,reffe_Disp;
          conformity=:H1,
          dirichlet_tags=["LeftSupp","LoadLine","RightSupp"],
          dirichlet_masks=[(true,true), (false,true),(false,true)])

uh = uh₈ = uh₉ = uh₁₀ = vh = ah = zero(V0_Disp)
initial_value = zero(V0_Disp)
	
end

# ╔═╡ 44072431-a452-4d11-a60c-3b79cceda615
function  stepPhaseField(uh_in,ψPlusPrev_in)
        
        a_PF(s,ϕ) = ∫( Gc*ls*∇(ϕ)⋅ ∇(s) + 2*ψPlusPrev_in*s*ϕ + (Gc/ls)*s*ϕ )*dΩ
        b_PF(ϕ) = ∫( (Gc/ls)*ϕ )*dΩ
        op_PF = AffineFEOperator(a_PF,b_PF,U_PF,V0_PF)
        sh_out = solve(op_PF)           
    
    return sh_out
    
end

# ╔═╡ c641c1d5-d941-41e4-b6a7-4c74f02eb2e0
function   stepDisp(uh_in, sh_in, count, time, T1, dt, v₀, uh, vh, ah)
        
        function uapp(x,t, time,T1, v₀)
           if time<=T1
              return VectorValue(0.0,-1.0)*((v₀/(2*T1))*t^2)
           else
              return VectorValue(0.0,-1.0)*((v₀*t) - (v₀*T1/2))
           end
        end
	
        u(x,t) = uapp(x,t, time,T1, v₀)
        g(x,t) = VectorValue(x[1]*0,x[2]*0)*t
        g(t::Real) = x-> g(x,t)
        u(t::Real) = x -> u(x,t)
    
        t0 = time-dt
        tF = time  
        Δt = dt/10
    
        U_Disp = TransientTrialFESpace(V0_Disp,[g,u,g])

        f(t) = VectorValue(0.0,0.0)
        m(utt,w) = ρ*∫(utt⋅w)dΩ
        c(ut,w) = ∫(ut⋅w)dΩ
            
        a_Disp(u,w) = ∫( (ε(w) ⊙ (σ_mod∘(ε(u),ε(uh_in),sh_in)) ) )*dΩ
        b_Disp(w,t) = ∫(w⋅f(t))dΩ
    
        res(t,u,w) = m(∂tt(u),w) + a_Disp(u,w) - b_Disp(w,t)
        jac(t,u,du,w) = a_Disp(du,w)
        jac_t(t,u,dut,w) = c(dut,w)
        jac_tt(t,u,dutt,w) = m(dutt,w)

        op = TransientFEOperator(res,jac,jac_t,jac_tt,U_Disp,V0_Disp)
    
        γ = 0.5
        β = 0.25

        U0 = U_Disp(t0)
        uh0 = interpolate_everywhere(uh,U0)
        vh0 = interpolate_everywhere(vh,U0)
        ah0 = interpolate_everywhere(ah,U0)

        Ls = LUSolver()
        ode_solver = Newmark(Ls,Δt,γ,β)

       sol_t = solve(ode_solver,op,(uh0,vh0,ah0),t0,tF) 
       
       counter = 0
        
       vApp =  u([L/2,H],time)[2]

       for (uₕ,t) in sol_t
          counter = counter + 1
          if counter == 8
            copyto!(uh₈.free_values, copy(uₕ.free_values))
            copyto!(uh₈.dirichlet_values, copy(uₕ.dirichlet_values))
          elseif counter == 9
            copyto!(uh₉.free_values, copy(uₕ.free_values)) 
            copyto!(uh₉.dirichlet_values, copy(uₕ.dirichlet_values))  
          elseif counter == 10
            copyto!(uh₁₀.free_values, copy(uₕ.free_values))
            copyto!(uh₁₀.dirichlet_values, copy(uₕ.dirichlet_values))
            print(" 👍 Doing at time $t") 
          end
       end  
      
       return (vApp, uh₈, uh₉, uh₁₀)
         
    end

# ╔═╡ a0c11cf1-10ce-4479-9be1-7f4edb9303c4
begin 
    sPrev = CellState(1.0,dΩ)
	sh = project(sPrev,model,dΩ,order)
	ψPlusPrev = CellState(0.0,dΩ)

	Load = Float64[]
    Displacement = Float64[]

    push!(Load, 0.0)
    push!(Displacement, 0.0)
end

# ╔═╡ 6e4b24d3-2eb7-4d3f-95c9-93ddac1280e8
md"""
## Defining the final dynamic-run simulation

"""

# ╔═╡ 2d7305c4-8503-4d44-a9d8-52a642cd0e1a
function dynamic_run(sh, uh, vh, ah)


tApp = 0
vApp = 0.0
v₀ = 60e-6
innerMax = 16
cnt = 0
dt1 =  4
Tmax = 2800
T1 = 196

	
while tApp .< Tmax 
    
    if tApp.<=T1 
        dt = dt1
    else
        dt = dt1*10
    end
    
    cnt = cnt .+ 1    
    tApp = tApp + dt
    
    print("\n Entering time step$cnt:", float(tApp))
    
    
   for inner = 1:innerMax
        
        ψhPlusPrev = project(ψPlusPrev,model,dΩ,order)
        
        RelErr = abs(sum(∫( Gc*ls*∇(sh)⋅ ∇(sh) + 2*ψhPlusPrev*sh*sh + (Gc/ls)*sh*sh)*dΩ-∫( (Gc/ls)*sh)*dΩ))/abs(sum(∫( (Gc/ls)*sh)*dΩ))
        print("\n Relative error = ",float(RelErr))
        
        sh = stepPhaseField(uh,ψhPlusPrev) 
        Δt = dt/10

        if cnt == 1
           vApp, uh8, uh9, uh = stepDisp(uh, sh, cnt, tApp, T1, dt, v₀, initial_value, initial_value, initial_value)
           vh9 = (uh9 - uh8)/Δt
           vh = (uh - uh9)/Δt
           ah = (vh - vh9)/Δt
        else
           vApp, uh8, uh9, uh = stepDisp(uh, sh, cnt, tApp, T1, dt, v₀, uh, vh, ah)
           vh9 = (uh9 - uh8)/Δt
           vh = (uh - uh9)/Δt
           ah = (vh - vh9)/Δt
        end

        ψhPos_in = ψPos∘(ε(uh))    
        update_state!(new_EnergyState,ψPlusPrev,ψhPos_in)
        
        if  RelErr < 1e-8
            break 
        end      
    end
    
    Node_Force = sum(∫( n_Γ_Load ⋅ (σ_mod∘(ε(uh),ε(uh),sh)) ) *dΓ_Load)
    push!(Load, -Node_Force[2])
    push!(Displacement, -vApp)
 
        writevtk(Ω,"results_JohnShah$cnt",cellfields=
         ["uh"=>uh,"s"=>sh , "epsi"=>ε(uh),"sigma"=>σ_mod∘(ε(uh),ε(uh),sh)])
end

end

# ╔═╡ bab63cba-55f1-42cc-b996-b1ad1ae5c3c7
dynamic_run(sh, uh, vh, ah)

# ╔═╡ 5263b040-3e32-45fb-b32c-23df937e2919
md"""
## Load v.s. Displacement curve

"""

# ╔═╡ e30b3aa5-f2c1-4281-a985-3be40ea0a572
plot(Displacement[1:end],Load[1:end].*1e-3, xlabel = "mm", ylabel = "kN")

# ╔═╡ Cell order:
# ╟─a2764d82-0497-4111-818e-9a0ffbb62796
# ╠═9892a47e-b7fb-11ed-1236-2f700b706c6b
# ╟─7a3de35c-6ba9-48ce-bbfb-e761a60cbbde
# ╠═601d363c-9fdd-4093-9fa5-b88ba7ce6cba
# ╟─35750bbc-d306-4806-9de0-0a3bbdaebcdb
# ╟─e36ffeea-13fa-431f-ae48-b44a6fedc76c
# ╠═b897e239-36ed-443a-beb8-965d46b1c0e5
# ╠═9a44e0aa-e92d-486d-920e-93918ae0f791
# ╠═08d10d55-bc28-4ac6-afb7-1d079a62548d
# ╟─13a77d13-7285-497f-a466-e0ac1685d9a5
# ╠═d39cc524-02b7-430d-ac2a-b0e2ffbe9e41
# ╠═af71ff02-1d8b-470a-bbf2-72174eeadeaa
# ╠═9621fd3c-15a3-4403-b805-f49245a16726
# ╟─1a22a2ba-3074-4b3a-9229-c0195645177a
# ╠═0af2040c-b6c4-42f3-848c-2def4b3e6881
# ╠═3d5dad6b-1db1-4dd6-b01f-23a5de21f9b0
# ╠═753a4974-059e-432f-b344-13504b9d401c
# ╠═1467f800-5ea9-4724-8749-81f2c4c467e7
# ╠═2a0d35c5-66e2-408c-b87d-f195eff673f6
# ╠═d8265d0e-a3e9-4a7f-8792-9275222b1810
# ╠═790de406-a97e-4988-896d-e8d62115cf36
# ╠═830e4d8e-8f22-4eed-9e0d-d0c3a943a72c
# ╟─b1d66ae6-5414-424e-b18d-c503082df7a8
# ╠═e8707054-6662-48a1-b758-3438990a5a48
# ╠═d2b2cb89-2d65-4a27-b34d-b8576218198e
# ╠═44072431-a452-4d11-a60c-3b79cceda615
# ╠═c641c1d5-d941-41e4-b6a7-4c74f02eb2e0
# ╠═a0c11cf1-10ce-4479-9be1-7f4edb9303c4
# ╟─6e4b24d3-2eb7-4d3f-95c9-93ddac1280e8
# ╠═2d7305c4-8503-4d44-a9d8-52a642cd0e1a
# ╠═bab63cba-55f1-42cc-b996-b1ad1ae5c3c7
# ╟─5263b040-3e32-45fb-b32c-23df937e2919
# ╠═e30b3aa5-f2c1-4281-a985-3be40ea0a572
