# ==============================================================================
# SANDER CMake File
# ==============================================================================

# program objects --------------------------------------------------------------
SET(SANDER_SRC
# -----------------------------------
        ../../lib/sebomd/se_deeprec.F90
        ../../lib/sebomd/se_info_from_sander.F90
        ../../lib/sebomd/se_mpi_vars.F90
        ../../lib/sebomd/se_inter.F90
        ../../lib/sebomd/se_corrections_info_from_sander.F90
        ../../lib/sebomd/se_corrections_module.F90
        ../../lib/sebomd/se_corrections_params.F90
        ../../lib/sebomd/se_corrections_tools.F90
        ../../lib/sebomd/se_write_charges.F90
        ../../lib/sebomd/se_edriver.F
        ../../lib/sebomd/se_corrections.F
        ../../lib/sebomd/sebomd.F
        ../../lib/sebomd/se_init.F90
        ../../lib/sebomd/se_energy.F
        ../../lib/sebomd/se_doscf.F
        ../../lib/sebomd/se_edriver.F
        ../../lib/sebomd/se_pme_setup.F
        ../../lib/sebomd/se_bond_order_wrapper.F
        ../../lib/sebomd/se_bond_order_write.F90

        ../../lib/sqm/qm2_dftb_scf.F90
        ../../lib/sqm/qm2_dftb_energy.F90 
        ../../lib/sqm/qm2_dftb_ewald.F90 
        ../../lib/sqm/qm2_diagonalizer_module.F90
        ../../lib/sqm/qm2_energy.F90
        ../../lib/sqm/qm2_get_qm_forces.F90 
        ../../lib/sqm/qm2_get_qmmm_forces.F90 
        ../../lib/sqm/qm2_load_params_and_allocate.F90 
        ../../lib/sqm/qm2_read_nm_and_alloc.F90
        ../../lib/sqm/qm2_scf.F90

# -----------------------------------
        abfqmmm_module.F90
        align.F90
#        ambmask.F90
        quench.F90
        relaxmd.F90
        amd.F90
        amoeba_adjust.F90
        amoeba_direct.F90
        amoeba_induced.F90
        amoeba_interface.F90
        amoeba_mdin.F90
        amoeba_multipoles.F90
        amoeba_recip.F90
        amoeba_runmd.F90
        amoeba_self.F90
        amoeba_valence.F90
        amoeba_vdw.F90
        amopen.F90
        ang_bend.F90
        apbs.F90
        apbs_vars.F90
        assert.F90
        barostats.F90
        binrestart.F90
        bintraj.F90
#       bmat_debug.F90
#        bmat_grad.F90
#        bmat_inv_debug.F90
#        bmat_inv.F90
#        bnd_strtch.F90
#        bond_anal2num.F90
#        cart2internal.F90
        charmm.F90
        check_input.F90
#        cmd_matrix_a1st.F90
        cmd_matrix.F90
        cmd_vars.F90
        constantph.F90
#        constants.F90
        crg_reloc.F90
        crossv.F90
        csa.F90
        cshf.F90
#        dbonds_anal2num.F90
        debug.F90
        decnvh.F90
        decomp.F90
        degcnt.F90
#        dg_debug.F90
#        dg_grad2cart.F90
#        dg_grad_anal2num.F90
        dihed_tors.F90
        dssp.F90
        dsvdc.F90
        dynlib.F90
#        egap_umb_2stdebug.F90
        egb.F90
        emap.F90
        emil.F90
        ene.F90
        erfcfun.F90
#        evb_2stdebug.F90
#        evb_alloc.F90
#        evb_bcast.F90
#        evb_dealloc.F90
#        evb_force.F90
#        evb_gfchk.F90
#        evb_init.F90
#        evb_input.F90
#        evb_io_debug2.F90
#        evb_io_debug.F90
#        evb_io.F90
#        evb_keywrd.F90
#        evb_mathf.F90
#        evb_matrix.F90
#        evb_ntrfc.F90
#        evb_outoptions.F90
#        evb_pimd_init.F90
#        evb_umb.F90
#        evb_umb_primitive.F90
#        evb_vars.F90
        ew_box.F90
        ew_bspline.F90
        ew_dipole_recip.F90
        ew_fft.F90
        ew_force.F90
        ew_handle_dips.F90
        ew_recip.F90
        ew_recip_reg.F90
        ew_setup.F90
#        exchange_gauss.F90
#        exchange_warshel.F90
        extra_pts.F90
        fastwt.F90
        file_io_dat.F90
        findmask.F90
#        fixport.F90
        force.F90
        getcor.F90
        icosasurf.F90
        interactmask.F90
        ips.F90
        iwrap2.F90
#        kappa_keck.F90
        linear_response.F90
#        lmod.F90
        les.F90
        lmod.F90
        locmem.F90
        lsc_init.F90
        lscivr_vars.F90
        lsc_xp.F90
        matinv.F90
        mbar.F90
        mdfil.F90
        mdin_emil_dat.F90
        mdread.F90
        mdwrit.F90
        memory_module.F90
        mexit.F90
        minrit.F90
#        mmtsb.F90
#        mod_vdw.F90
        molecule.F90
#        morse_anal2num.F90
#        morsify.F90
        mt19937.F90
        mtmdcall.F90
        multisander.F90
        multitmd.F90
nfe-abmd-hooks.F90
nfe-bbmd-ctxt.F90
nfe-bbmd-hooks.F90
nfe-colvar.F90
nfe-colvar-math.F90
nfe-colvar-type.F90
nfe-colvar-utils.F90
nfe-constants.F90
nfe-cv-ANGLE.F90
nfe-cv-COM_ANGLE.F90
nfe-cv-COM_DISTANCE.F90
nfe-cv-COM_TORSION.F90
nfe-cv-COS_OF_DIHEDRAL.F90
nfe-cv-DF_COM_DISTANCE.F90
nfe-cv-DISTANCE.F90
nfe-cv-HANDEDNESS.F90
nfe-cv-LCOD.F90
nfe-cv-MULTI_RMSD.F90
nfe-cv-N_OF_BONDS.F90
nfe-cv-N_OF_STRUCTURES.F90
nfe-cv-PAIR_DIHEDRAL.F90
nfe-cv-PATTERN_DIHEDRAL.F90
nfe-cv-PCA.F90
nfe-cv-R_OF_GYRATION.F90
nfe-cv-SIN_OF_DIHEDRAL.F90
nfe-cv-TORSION.F90
nfe-pmd-hooks.F90
nfe-read-pca.F90
nfe-rmsd.F90
nfe-sander-hooks.F90
nfe-sander-proxy.F90
nfe-smd-hooks.F90
nfe-stsm-hooks.F90
nfe-umbrella.F90
nfe-utils.F90
nfe-value.F90
        new_time.F90
        nmlsrc.F90
        nmrcal.F90
        nmr.F90
        nonbond_list.F90
        nose_hoover.F90
        nose_hoover_init.F90
        nose_hoover_init_LES.F90
#        out_evb.F90
        parallel.F90
        parms.F90
        pcshift.F90
        pearsn.F90
        pimd_force.F90
        pimd_init.F90
        pimd_vars.F90
        printe.F90
        prn_dipoles.F90
        prn_qmmm_dipole.F90
#        pubfft.F90
        pupildata.F90
        putvalues.F90
        qm2_extern_adf_module.F90
        qm2_read_adf_results.c
        KFReader.c
        ArrayList.c
        qm2_extern_gau_module.F90
        qm2_extern_genmpi_module.F90
        qm2_extern_gms_module.F90
#        qm2_extern_lio_module.F90
        qm2_extern_module.F90
        qm2_extern_nw_module.F90
        qm2_extern_orc_module.F90
        qm2_extern_tc_module.F90
        qm2_extern_util_module.F90
        qm2_extern_qc_module.F90
        qm2_extern_mrcc_module.F90
#        qm_div.F90
        qm_ewald.F90
        qmmm_adaptive_module.F90
        qm_mm.F90
        qmmm_vsolv.F90
#        rcbond_anal2num.F90
#        rcdbonds_anal2num.F90
        rdparm.F90
#        react_flux.F90
#        react_flux_init.F90
#        react_path.F90
        relax_mat.F90
        remd.F90
        rfree.F90
        rgroup.F90
        rmsgrd.F90
        runmd.F90
        runmin.F90
#        runext.F90
        sander.F90
        sander_lib.F90
        scaledMD.F90
#        schlegel_angle.F90
#        schlegel_dg.F90
#        schlegel_dihed.F90
#        schlegel_diis.F90
#        schlegel_full.F90
#        schlegel_gbasis.F90
#        schlegel_gmres.F90
#        schlegel_hcn_cart.F90
#        schlegel_hcn_cart_irc.F90
#        schlegel_hcn_int_irc.F90
#        schlegel_pgbasis.F90
#        schlegel_poh_irc.F90
#        schlegel_poh_uff.F90
#        schlegel_uff.F90
        sebomd_arrays.F90
        sebomd_module.F90
#        select_coord.F90
        set.F90
        sglds.F90
        shake.F90
        short_ene.F90
        softcore.F90
        spatial_fft.F90
        spatial_recip.F90
        stack.F90
        state.F90
        thermo_int.F90
#        torinternal.F90
        trace.F90
        trajene.F90
        unitv.F90
#        vdw_anal2num.F90
#        wdc_bmat.F90
#        xgauss_anal2num.F90
        xray_fftpack.F90
        xray_fourier.F90
        xray_globals.F90
        xray_interface.F90
        xray_real_space.F90
        xray_reciprocal_space.F90
        xray_utils.F90
        xref.F90
#        xwarshel_anal2num.F90
#        yammpnb.F90
        Plumed.c
        pol_gauss.F90
        ew_bndpol.F90
        sinr_t.F90
        )


IF(MPI_ENABLED)
    SET(SANDER_SRC ${SANDER_SRC}
    # -----------------------------------
        bmat_debug.F90
        bmat_grad.F90
        bmat_inv_debug.F90
        bmat_inv.F90
        bnd_strtch.F90
        bond_anal2num.F90
        cart2internal.F90
#        cmd_matrix_a1st.F90
        dbonds_anal2num.F90
        dg_debug.F90
        dg_grad2cart.F90
        dg_grad_anal2num.F90
        egap_umb_2stdebug.F90
        evb_2stdebug.F90
        evb_alloc.F90
        evb_bcast.F90
        evb_dealloc.F90
        evb_force.F90
        evb_gfchk.F90
        evb_init.F90
        evb_input.F90
        evb_io_debug2.F90
        evb_io_debug.F90
        evb_io.F90
        evb_keywrd.F90
        evb_mathf.F90
        evb_matrix.F90
        evb_ntrfc.F90
        evb_outoptions.F90
        evb_pimd_init.F90
        evb_umb.F90
        evb_umb_primitive.F90
        evb_vars.F90
        exchange_gauss.F90
        exchange_warshel.F90
#        fixport.F90
        kappa_keck.F90
        mod_vdw.F90
        morse_anal2num.F90
        morsify.F90
        out_evb.F90
#        pubfft.F90
        qm_div.F90
        rcbond_anal2num.F90
        rcdbonds_anal2num.F90
        react_flux.F90
        react_flux_init.F90
        react_path.F90
        schlegel_angle.F90
        schlegel_dg.F90
        schlegel_dihed.F90
        schlegel_diis.F90
        schlegel_full.F90
        schlegel_gbasis.F90
        schlegel_gmres.F90
        schlegel_hcn_cart.F90
        schlegel_hcn_cart_irc.F90
        schlegel_hcn_int_irc.F90
        schlegel_pgbasis.F90
        schlegel_poh_irc.F90
        schlegel_poh_uff.F90
        schlegel_uff.F90
        select_coord.F90
        torinternal.F90
        vdw_anal2num.F90
        wdc_bmat.F90
        xgauss_anal2num.F90
        xwarshel_anal2num.F90
        )
ENDIF(MPI_ENABLED)

# ------------------------------------------------------------------------------

ADD_DEFINITIONS(${SANDER_DEFS})

IF(MPI_ENABLED)
    SET(TARGET_NAME _sander.PMF.MPI)
ELSE(MPI_ENABLED)
    SET(TARGET_NAME sander.PMF)
ENDIF(MPI_ENABLED)

ADD_EXECUTABLE(${TARGET_NAME} ${SANDER_SRC})

IF(HAVE_FFTW3)
	ADD_DEPENDENCIES(${TARGET_NAME} utils_shared pbsa_shared sqm_shared sebomd_shared rism_shared)
ELSE(HAVE_FFTW3)
	ADD_DEPENDENCIES(${TARGET_NAME} utils_shared pbsa_shared sqm_shared sebomd_shared)
ENDIF(HAVE_FFTW3)	

TARGET_LINK_LIBRARIES(${TARGET_NAME}
        ${SANDER_LIBRARIES}
        ${FPMF_LIBRARIES}
        )

INSTALL(TARGETS
            ${TARGET_NAME}
        DESTINATION
            bin
        )

# ------------------------------------------------------------------------------
