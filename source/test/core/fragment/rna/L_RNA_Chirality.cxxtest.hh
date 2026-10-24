// -*- mode:c++;tab-width:2;indent-tabs-mode:t;show-trailing-whitespace:t;rm-trailing-spaces:t -*-
// vi: set ts=2 noet:
//
// (c) Copyright Rosetta Commons Member Institutions.
// (c) This file is part of the Rosetta software suite and is made available under license.
// (c) The Rosetta software is developed by the contributing members of the Rosetta Commons.
// (c) For more information, see http://www.rosettacommons.org. Questions about this can be
// (c) addressed to University of Washington CoMotion, email: license@uw.edu.

/// @file   test/core/fragment/rna/L_RNA_Chirality.cxxtest.hh
/// @brief  Mirror-image (L-RNA) handling: fragment insertion, pucker helpers, ideal geometry, scoring.
/// @details Each test builds a D-RNA pose and its exact mirror image out of the built-in L-RNA
///  types (one-letter codes 1-4), then checks that L-RNA code paths treat the mirror image exactly
///  as the D-RNA paths treat the original.

// Test headers
#include <cxxtest/TestSuite.h>
#include <test/core/init_util.hh>

#include <core/types.hh>
#include <core/chemical/AA.hh>
#include <core/chemical/ResidueType.hh>
#include <core/chemical/rna/util.hh>
#include <core/conformation/Conformation.hh>
#include <core/conformation/Residue.hh>
#include <core/fragment/rna/FullAtomRNA_Fragments.hh>
#include <core/fragment/rna/TorsionSet.hh>
#include <core/id/AtomID.hh>
#include <core/id/TorsionID.hh>
#include <core/pose/Pose.hh>
#include <core/pose/annotated_sequence.hh>
#include <core/pose/rna/util.hh>
#include <core/pose/rna/RNA_IdealCoord.hh>
#include <core/pose/toolbox/AtomLevelDomainMap.hh>
#include <core/scoring/Energies.hh>
#include <core/scoring/ScoreFunction.hh>
#include <core/scoring/ScoreFunctionFactory.hh>
#include <core/scoring/ScoreType.hh>
#include <core/sequence/util.hh>
#include <core/chemical/ResidueTypeSet.hh>
#include <core/chemical/VariantType.hh>
#include <core/import_pose/libraries/BasePairStepLibrary.hh>
#include <core/kinematics/FoldTree.hh>
#include <core/kinematics/MoveMap.hh>
#include <core/pose/MiniPose.hh>
#include <core/pose/rna/BasePairStep.hh>
#include <core/pose/variant_util.hh>
#include <core/pose/full_model_info/FullModelInfo.hh>
#include <core/chemical/ChemicalManager.hh>
#include <core/optimization/AtomTreeMultifunc.hh>
#include <core/optimization/MinimizerMap.hh>
#include <core/optimization/NumericalDerivCheckResult.hh>
#include <core/optimization/types.hh>
#include <core/optimization/atom_tree_minimize.hh>

#include <basic/database/open.hh>
#include <basic/Tracer.hh>

#include <cmath>

static basic::Tracer TR( "core.fragment.rna.L_RNA_Chirality.cxxtest" );

using namespace core;
using namespace core::chemical::rna;

class L_RNA_ChiralityTests : public CxxTest::TestSuite {

public:

	void setUp() {
		core_init();
	}

	void tearDown() {}

	// All RNA torsions (backbone + chi, nu2, nu1, 2'-OH) that exist in this pose.
	utility::vector1< id::TorsionID >
	rna_torsion_ids( pose::Pose const & pose ) {
		utility::vector1< id::TorsionID > torsion_ids;
		for ( Size i = 1; i <= pose.size(); ++i ) {
			for ( Size j = 1; j <= NUM_RNA_TORSIONS; ++j ) {
				id::TorsionID const tid = ( j <= NUM_RNA_MAINCHAIN_TORSIONS ) ?
					id::TorsionID( i, id::BB, j ) : id::TorsionID( i, id::CHI, j - NUM_RNA_MAINCHAIN_TORSIONS );
				id::AtomID id1, id2, id3, id4;
				if ( pose.conformation().get_torsion_angle_atom_ids( tid, id1, id2, id3, id4 ) ) continue; // true means failure
				torsion_ids.push_back( tid );
			}
		}
		return torsion_ids;
	}

	// Make l_pose the mirror image of d_pose by negating every torsion.
	void
	mirror_torsions( pose::Pose const & d_pose, pose::Pose & l_pose ) {
		for ( id::TorsionID const & tid : rna_torsion_ids( d_pose ) ) {
			l_pose.set_torsion( tid, -d_pose.torsion( tid ) );
		}
	}

	// Mirror images have identical internal distances.
	void
	assert_mirror_images( pose::Pose const & d_pose, pose::Pose const & l_pose, Real const tol = 1.0e-3 ) {
		TS_ASSERT_EQUALS( d_pose.size(), l_pose.size() );
		for ( Size i = 1; i <= d_pose.size(); ++i ) {
			TS_ASSERT_EQUALS( d_pose.residue( i ).natoms(), l_pose.residue( i ).natoms() );
			for ( Size j = i; j <= d_pose.size(); ++j ) {
				for ( Size ii = 1; ii <= d_pose.residue( i ).nheavyatoms(); ++ii ) {
					for ( Size jj = 1; jj <= d_pose.residue( j ).nheavyatoms(); ++jj ) {
						Real const d_dist = d_pose.residue( i ).xyz( ii ).distance( d_pose.residue( j ).xyz( jj ) );
						Real const l_dist = l_pose.residue( i ).xyz( ii ).distance( l_pose.residue( j ).xyz( jj ) );
						TS_ASSERT_DELTA( d_dist, l_dist, tol );
					}
				}
			}
		}
	}

	void test_one_letter_codes_give_l_rna() {
		pose::Pose l_pose;
		pose::make_pose_from_sequence( l_pose, "1234", "fa_standard" );
		// L-RNA types keep the RNA letters as one-letter codes; 1-4 is input shorthand only.
		TS_ASSERT_EQUALS( l_pose.sequence(), "acgu" );
		for ( Size i = 1; i <= l_pose.size(); ++i ) TS_ASSERT( l_pose.residue_type( i ).is_l_rna() );
		TS_ASSERT_EQUALS( l_pose.residue_type( 1 ).aa(), chemical::na_lra );
		TS_ASSERT_EQUALS( l_pose.residue_type( 2 ).aa(), chemical::na_lrc );
		TS_ASSERT_EQUALS( l_pose.residue_type( 3 ).aa(), chemical::na_lrg );
		TS_ASSERT_EQUALS( l_pose.residue_type( 4 ).aa(), chemical::na_lur );

		// FASTA / -sequence input: shorthand expands to the annotated form.
		TS_ASSERT_EQUALS( sequence::expand_l_rna_shorthand( "ga13u" ), "gaa[LRAD]g[LRGU]u" );
		TS_ASSERT_EQUALS( sequence::expand_l_rna_shorthand( "g1u" ), "ga[LRAD]u" );
		TS_ASSERT_EQUALS( sequence::expand_l_rna_shorthand( "g3[LRGU:LowerRNA]u" ), "gg[LRGU:LowerRNA]u" );
		TS_ASSERT_EQUALS( sequence::expand_l_rna_shorthand( "Z[MG]a2" ), "Z[MG]ac[LRCY]" );
		std::string seq( "gc12" );
		std::map< Size, std::string > const nonstandard = sequence::parse_out_non_standard_residues( seq );
		TS_ASSERT_EQUALS( seq, "gcac" );
		TS_ASSERT_EQUALS( nonstandard.size(), 2 );
		TS_ASSERT_EQUALS( nonstandard.at( 3 ), "LRAD" );
		TS_ASSERT_EQUALS( nonstandard.at( 4 ), "LRCY" );
	}

	void test_mirrored_types_and_torsions_give_mirror_image() {
		pose::Pose d_pose, l_pose;
		pose::make_pose_from_sequence( d_pose, "gcaugc", "fa_standard" );
		pose::make_pose_from_sequence( l_pose, "321432", "fa_standard" );
		mirror_torsions( d_pose, l_pose );
		assert_mirror_images( d_pose, l_pose );
	}

	void test_pucker_assignment_and_ideal_c2endo() {
		pose::Pose d_pose, l_pose;
		pose::make_pose_from_sequence( d_pose, "gcaugc", "fa_standard" );
		pose::make_pose_from_sequence( l_pose, "321432", "fa_standard" );
		pose::rna::apply_ideal_c2endo_sugar_coords( d_pose, 3 );
		pose::rna::apply_ideal_c2endo_sugar_coords( l_pose, 3 );
		TS_ASSERT_DELTA( l_pose.torsion( id::TorsionID( 3, id::BB, DELTA ) ), -d_pose.torsion( id::TorsionID( 3, id::BB, DELTA ) ), 1.0e-3 );
		for ( Size i = 1; i <= d_pose.size(); ++i ) {
			TS_ASSERT_EQUALS( pose::rna::assign_pucker( d_pose, i ), pose::rna::assign_pucker( l_pose, i ) );
			TS_ASSERT_EQUALS( get_residue_pucker_state( d_pose.residue( i ) ), get_residue_pucker_state( l_pose.residue( i ) ) );
			TS_ASSERT_EQUALS( get_residue_base_state( d_pose.residue( i ) ), get_residue_base_state( l_pose.residue( i ) ) );
		}
		TS_ASSERT_EQUALS( pose::rna::assign_pucker( l_pose, 3 ), SOUTH );
		// The rest of each pose is still the default A-form, so the two are full mirror images.
		mirror_torsions( d_pose, l_pose );
		assert_mirror_images( d_pose, l_pose );
	}

	void test_ideal_coord_pucker_is_mirrored() {
		pose::Pose d_pose, l_pose;
		pose::make_pose_from_sequence( d_pose, "gcaugc", "fa_standard" );
		pose::make_pose_from_sequence( l_pose, "321432", "fa_standard" );
		mirror_torsions( d_pose, l_pose );
		pose::rna::RNA_IdealCoord const ideal_coord;
		for ( Size i = 1; i <= d_pose.size(); ++i ) {
			ideal_coord.apply( d_pose, i, SOUTH );
			ideal_coord.apply( l_pose, i, SOUTH );
		}
		for ( Size i = 1; i <= d_pose.size(); ++i ) {
			TS_ASSERT_EQUALS( pose::rna::assign_pucker( d_pose, i ), SOUTH );
			TS_ASSERT_EQUALS( pose::rna::assign_pucker( l_pose, i ), SOUTH );
		}
		assert_mirror_images( d_pose, l_pose );
	}

	void test_fragment_insertion_inverts_for_l_and_skips_mixed_windows() {
		using namespace core::fragment::rna;
		FullAtomRNA_Fragments const fragments( basic::database::full_name( "sampling/rna/1jj2.torsions" ) );

		// A two-residue fragment with a C2'-endo first residue and non-A-form backbone.
		TorsionSet torsion_set( 2, 1 );
		Real const values[ 2 ][ NUM_RNA_TORSIONS ] = {
			{ 69.4, -173.0, 58.9, 147.2, -85.4, -38.4, 111.7, -36.4, 156.4, 179.9 },
			{ -68.0, 178.0, 54.0, 82.0, -150.0, -71.0, 79.0, 37.0, -25.0, 100.0 } };
		for ( Size offset = 0; offset < 2; ++offset ) {
			for ( Size j = 1; j <= NUM_RNA_TORSIONS; ++j ) torsion_set.torsions( j, offset ) = values[ offset ][ j - 1 ];
		}

		pose::Pose d_pose, l_pose, mixed_pose;
		pose::make_pose_from_sequence( d_pose, "gcaugc", "fa_standard" );
		pose::make_pose_from_sequence( l_pose, "321432", "fa_standard" );
		pose::make_pose_from_sequence( mixed_pose, "gca432", "fa_standard" );

		pose::toolbox::AtomLevelDomainMapOP d_map( new pose::toolbox::AtomLevelDomainMap( d_pose ) );
		pose::toolbox::AtomLevelDomainMapOP l_map( new pose::toolbox::AtomLevelDomainMap( l_pose ) );
		pose::toolbox::AtomLevelDomainMapOP mixed_map( new pose::toolbox::AtomLevelDomainMap( mixed_pose ) );
		d_map->set( true ); l_map->set( true ); mixed_map->set( true );

		bool is_l( false );
		TS_ASSERT( window_has_uniform_rna_chirality( d_pose, 3, 2, is_l ) && !is_l );
		TS_ASSERT( window_has_uniform_rna_chirality( l_pose, 3, 2, is_l ) && is_l );
		TS_ASSERT( !window_has_uniform_rna_chirality( mixed_pose, 3, 2, is_l ) );

		fragments.insert_fragment( d_pose, 3, torsion_set, d_map );
		fragments.insert_fragment( l_pose, 3, torsion_set, l_map );
		for ( Size offset = 0; offset < 2; ++offset ) {
			for ( Size j = 1; j <= NUM_RNA_MAINCHAIN_TORSIONS; ++j ) {
				id::TorsionID const tid( 3 + offset, id::BB, j );
				TS_ASSERT_DELTA( d_pose.torsion( tid ), values[ offset ][ j - 1 ], 1.0e-3 );
				TS_ASSERT_DELTA( l_pose.torsion( tid ), -values[ offset ][ j - 1 ], 1.0e-3 );
			}
			id::TorsionID const chi( 3 + offset, id::CHI, 1 );
			TS_ASSERT_DELTA( l_pose.torsion( chi ), -d_pose.torsion( chi ), 1.0e-3 );
		}
		mirror_torsions( d_pose, l_pose ); // only touches torsions; should already hold
		assert_mirror_images( d_pose, l_pose );

		// Mixed window: nothing changes.
		pose::Pose const mixed_before( mixed_pose );
		fragments.insert_fragment( mixed_pose, 3, torsion_set, mixed_map );
		for ( id::TorsionID const & tid : rna_torsion_ids( mixed_pose ) ) {
			TS_ASSERT_DELTA( mixed_pose.torsion( tid ), mixed_before.torsion( tid ), 1.0e-6 );
		}
	}

	// Scoring a mirror-image L-RNA pose must give the same energies, term by term, as the D pose.
	void
	check_mirror_energies( std::string const & weights ) {
		scoring::ScoreFunctionOP sfxn = scoring::ScoreFunctionFactory::create_score_function( weights );
		sfxn->set_weight( scoring::rna_suite, 1.0 );
		sfxn->set_weight( scoring::rna_jr_suite, 1.0 );
		sfxn->set_weight( scoring::rna_sugar_close, 1.0 );
		sfxn->set_weight( scoring::rna_torsion, 1.0 );

		pose::Pose d_pose, l_pose;
		pose::make_pose_from_sequence( d_pose, "ggcuucggcc", "fa_standard" );
		pose::make_pose_from_sequence( l_pose, "3324423322", "fa_standard" );
		// Some non-A-form structure so that suites, puckers and stacking are not all identical.
		// (Applied to both poses: the C2'-endo sugar ring has non-ideal bond geometry, not just torsions.)
		for ( Size const i : { 5, 6 } ) {
			pose::rna::apply_ideal_c2endo_sugar_coords( d_pose, i );
			pose::rna::apply_ideal_c2endo_sugar_coords( l_pose, i );
		}
		d_pose.set_torsion( id::TorsionID( 6, id::BB, ALPHA ), 70.0 );
		d_pose.set_torsion( id::TorsionID( 7, id::CHI, 1 ), -60.0 ); // syn
		mirror_torsions( d_pose, l_pose );
		assert_mirror_images( d_pose, l_pose );

		Real const d_total = ( *sfxn )( d_pose );
		Real const l_total = ( *sfxn )( l_pose );
		for ( scoring::ScoreType const st : sfxn->get_nonzero_weighted_scoretypes() ) {
			Real const d_e = d_pose.energies().total_energies()[ st ];
			Real const l_e = l_pose.energies().total_energies()[ st ];
			TR << weights << " " << scoring::name_from_score_type( st ) << " D: " << d_e << " L: " << l_e << std::endl;
			TS_ASSERT_DELTA( l_e, d_e, 1.0e-3 * std::max( 1.0, std::abs( d_e ) ) );
		}
		TS_ASSERT_DELTA( l_total, d_total, 1.0e-3 * std::max( 1.0, std::abs( d_total ) ) );
	}

	void test_mirror_energies_hires() {
		check_mirror_energies( "stepwise/rna/rna_res_level_energy4" );
	}

	void test_mirror_energies_lores() {
		check_mirror_energies( "rna/denovo/rna_lores" );
	}


	// A D pose and its exact L mirror, with some non-A-form geometry.
	void
	make_mirror_pair( pose::Pose & d_pose, pose::Pose & l_pose ) {
		pose::make_pose_from_sequence( d_pose, "ggcuucggcc", "fa_standard" );
		pose::make_pose_from_sequence( l_pose, "3324423322", "fa_standard" );
		for ( Size const i : { 5, 6 } ) {
			pose::rna::apply_ideal_c2endo_sugar_coords( d_pose, i );
			pose::rna::apply_ideal_c2endo_sugar_coords( l_pose, i );
		}
		d_pose.set_torsion( id::TorsionID( 6, id::BB, ALPHA ), 70.0 );
		d_pose.set_torsion( id::TorsionID( 7, id::CHI, 1 ), -60.0 ); // syn
		mirror_torsions( d_pose, l_pose );
	}

	// Numeric and analytic rna_torsion derivatives for every backbone/chi DOF of a pose (same DOF order for poses
	// with the same atom tree). Upstream rna_torsion has no analytic derivative for the zeta DOF (O5' phi) even in
	// D-RNA, so the tests below compare L with D rather than requiring numeric == analytic everywhere.
	utility::vector1< std::pair< Real, Real > >
	rna_torsion_derivs( pose::Pose pose ) {
		using namespace core::optimization;
		scoring::ScoreFunction sfxn;
		sfxn.set_weight( scoring::rna_torsion, 1.0 );
		kinematics::MoveMap mm;
		mm.set_bb( true );
		mm.set_chi( true );
		sfxn( pose );
		MinimizerMap min_map;
		min_map.setup( pose, mm );
		pose.energies().set_use_nblist( pose, min_map.domain_map(), false );
		sfxn.setup_for_minimizing( pose, min_map );
		sfxn.setup_for_derivatives( pose );
		AtomTreeMultifunc f( pose, min_map, sfxn, false, false );
		Multivec vars( min_map.nangles() );
		min_map.copy_dofs_from_pose( pose, vars );
		Multivec dE_dvars( vars );
		f.dfunc( vars, dE_dvars );
		SimpleDerivCheckResult const result = simple_numeric_deriv_check( f, vars, dE_dvars, false, false, 1 );
		utility::vector1< std::pair< Real, Real > > derivs;
		for ( int ii = 1; ii <= min_map.nangles(); ++ii ) {
			derivs.push_back( std::make_pair( result.step_data( ii, 1 ).num_deriv(), result.step_data( ii, 1 ).ana_deriv() ) );
		}
		return derivs;
	}

	// L-RNA derivatives are the mirror of D (the L f2 sign used to be wrong): analytic dE/dtorsion negated,
	// and exactly the same numeric-vs-analytic agreement as D.
	void test_rna_torsion_derivatives_d_and_l() {
		pose::Pose d_pose, l_pose;
		make_mirror_pair( d_pose, l_pose );
		utility::vector1< std::pair< Real, Real > > const d = rna_torsion_derivs( d_pose );
		utility::vector1< std::pair< Real, Real > > const l = rna_torsion_derivs( l_pose );
		TS_ASSERT_EQUALS( d.size(), l.size() );
		Size n_nonzero( 0 );
		for ( Size ii = 1; ii <= std::min( d.size(), l.size() ); ++ii ) {
			TS_ASSERT_DELTA( l[ ii ].second, -d[ ii ].second, 1.0e-6 );
			TS_ASSERT_DELTA( l[ ii ].first - l[ ii ].second, -( d[ ii ].first - d[ ii ].second ), 1.0e-4 );
			if ( std::abs( d[ ii ].second ) > 1.0e-3 ) ++n_nonzero;
		}
		TS_ASSERT( n_nonzero > 10 ); // the check is not vacuous
	}

	// A covalent D/L junction scores 0 in rna_torsion; its derivative terms must vanish too, so the junction adds
	// no numeric-vs-analytic mismatch beyond what an all-D chain with the same geometry has.
	void test_rna_torsion_derivatives_mixed_junction() {
		pose::Pose mixed_pose, d_pose;
		pose::make_pose_from_sequence( mixed_pose, "ggcc3322", "fa_standard" );
		pose::make_pose_from_sequence( d_pose, "ggccggcc", "fa_standard" );
		for ( Size i = 5; i <= 8; ++i ) {
			for ( Size j = 1; j <= 6; ++j ) d_pose.set_torsion( id::TorsionID( i, id::BB, j ), -mixed_pose.torsion( id::TorsionID( i, id::BB, j ) ) );
			d_pose.set_torsion( id::TorsionID( i, id::CHI, 1 ), -mixed_pose.torsion( id::TorsionID( i, id::CHI, 1 ) ) );
		}
		utility::vector1< std::pair< Real, Real > > const m = rna_torsion_derivs( mixed_pose );
		utility::vector1< std::pair< Real, Real > > const d = rna_torsion_derivs( d_pose );
		TS_ASSERT_EQUALS( m.size(), d.size() );
		for ( Size ii = 1; ii <= std::min( m.size(), d.size() ); ++ii ) {
			TS_ASSERT_LESS_THAN_EQUALS( std::abs( m[ ii ].first - m[ ii ].second ), std::abs( d[ ii ].first - d[ ii ].second ) + 1.0e-4 );
		}
	}

	// Base-pair-step chunks: the mirrored library is the exact reflection; step chirality is classified.
	void test_base_pair_step_mirroring() {
		using namespace core::import_pose::libraries;
		BasePairStepLibrary const library( true );
		BasePairStepSequence const seq( 'g', 'c', 'g', 'c' );
		TS_ASSERT( library.has_value( seq ) );
		utility::vector1< pose::MiniPoseOP > const & d_list = library.mini_pose_list( seq );
		utility::vector1< pose::MiniPoseOP > const & l_list = library.mirrored_mini_pose_list( seq );
		TS_ASSERT( d_list.size() > 0 );
		TS_ASSERT_EQUALS( d_list.size(), l_list.size() );
		for ( Size k = 1; k <= std::min( d_list.size(), Size( 3 ) ); ++k ) {
			for ( Size i = 1; i <= d_list[ k ]->coords().size(); ++i ) {
				for ( Size a = 1; a <= d_list[ k ]->coords()[ i ].size(); ++a ) {
					PointPosition const & d = d_list[ k ]->coords()[ i ][ a ];
					PointPosition const & l = l_list[ k ]->coords()[ i ][ a ];
					TS_ASSERT_DELTA( l.x(), -d.x(), 1.0e-9 );
					TS_ASSERT_DELTA( l.y(), d.y(), 1.0e-9 );
					TS_ASSERT_DELTA( l.z(), d.z(), 1.0e-9 );
				}
			}
		}
		pose::Pose pose;
		pose::make_pose_from_sequence( pose, "gggg3333gggg", "fa_standard" );
		bool is_l( false );
		TS_ASSERT( pose::rna::base_pair_step_has_uniform_rna_chirality( pose::rna::BasePairStep( 1, 2, 11, 12 ), pose, is_l ) );
		TS_ASSERT( !is_l );
		TS_ASSERT( pose::rna::base_pair_step_has_uniform_rna_chirality( pose::rna::BasePairStep( 5, 6, 7, 8 ), pose, is_l ) );
		TS_ASSERT( is_l );
		TS_ASSERT( !pose::rna::base_pair_step_has_uniform_rna_chirality( pose::rna::BasePairStep( 4, 5, 11, 12 ), pose, is_l ) );
	}

	// Canonical-pair identity: L:L pairs are canonical like D:D, D:L pairs are not.
	void test_canonical_pair_identity() {
		using namespace core::chemical;
		using namespace core::chemical::rna;
		TS_ASSERT( possibly_canonical( na_lrg, na_lrc ) );
		TS_ASSERT( possibly_canonical( na_lra, na_lur ) );
		TS_ASSERT( possibly_canonical( na_lrg, na_lur ) );
		TS_ASSERT( !possibly_canonical( na_rgu, na_lrc ) );
		TS_ASSERT( possibly_canonical_strict( na_lrg, na_lrc ) );
		TS_ASSERT( !possibly_canonical_strict( na_lrg, na_lur ) );
		TS_ASSERT_EQUALS( rna_d_equivalent_aa( na_lra ), na_rad );
		TS_ASSERT_EQUALS( rna_d_equivalent_aa( na_rgu ), na_rgu );

		ResidueTypeSetCOP rts = ChemicalManager::get_instance()->residue_type_set( FA_STANDARD );
		utility::vector1< std::string > d1, d2, l1, l2, m1, m2;
		get_base_pair_atoms( rts->name_map( "RGU" ), rts->name_map( "RCY" ), d1, d2, WATSON_CRICK, WATSON_CRICK, CIS );
		get_base_pair_atoms( rts->name_map( "LRGU" ), rts->name_map( "LRCY" ), l1, l2, WATSON_CRICK, WATSON_CRICK, CIS );
		get_base_pair_atoms( rts->name_map( "RGU" ), rts->name_map( "LRCY" ), m1, m2, WATSON_CRICK, WATSON_CRICK, CIS );
		TS_ASSERT( d1.size() > 0 );
		TS_ASSERT_EQUALS( d1, l1 );
		TS_ASSERT_EQUALS( d2, l2 );
		TS_ASSERT_EQUALS( m1.size(), 0 );
	}

	// Cutpoint variants (connection ICOORs, phosphate placement) give mirror-image virtual atoms.
	void test_cutpoint_variants_are_mirrored() {
		pose::Pose d_pose, l_pose;
		make_mirror_pair( d_pose, l_pose );
		for ( pose::Pose * pose : { &d_pose, &l_pose } ) {
			kinematics::FoldTree ft( pose->size() );
			ft.new_jump( 2, 9, 5 );
			pose->fold_tree( ft );
			pose::correctly_add_cutpoint_variants( *pose, 5 );
		}
		assert_mirror_images( d_pose, l_pose );
		for ( std::string const name : { "OVL1", "OVL2" } ) {
			for ( std::string const partner : { "P", "O5'" } ) {
				TS_ASSERT_DELTA( d_pose.residue( 5 ).xyz( name ).distance( d_pose.residue( 6 ).xyz( partner ) ),
					l_pose.residue( 5 ).xyz( name ).distance( l_pose.residue( 6 ).xyz( partner ) ), 1.0e-3 );
			}
		}
		TS_ASSERT_DELTA( d_pose.residue( 6 ).xyz( "OVU1" ).distance( d_pose.residue( 5 ).xyz( "O3'" ) ),
			l_pose.residue( 6 ).xyz( "OVU1" ).distance( l_pose.residue( 5 ).xyz( "O3'" ) ), 1.0e-3 );
		TS_ASSERT_DELTA( d_pose.residue( 6 ).xyz( "OP1" ).distance( d_pose.residue( 5 ).xyz( "C3'" ) ),
			l_pose.residue( 6 ).xyz( "OP1" ).distance( l_pose.residue( 5 ).xyz( "C3'" ) ), 1.0e-3 );
	}

	// L-RNA given by name at chain ends gets terminus variants, like one-letter D-RNA.
	void test_named_l_rna_gets_terminus_variants() {
		pose::Pose d_pose, l_pose;
		pose::make_pose_from_sequence( d_pose, "ggcc", "fa_standard" );
		pose::make_pose_from_sequence( l_pose, "g[LRGU]g[LRGU]c[LRCY]c[LRCY]", "fa_standard" );
		TS_ASSERT( d_pose.residue_type( 1 ).has_variant_type( chemical::LOWER_TERMINUS_VARIANT ) );
		TS_ASSERT( l_pose.residue_type( 1 ).has_variant_type( chemical::LOWER_TERMINUS_VARIANT ) );
		TS_ASSERT( l_pose.residue_type( 4 ).has_variant_type( chemical::UPPER_TERMINUS_VARIANT ) );
		TS_ASSERT( !l_pose.residue_type( 2 ).has_variant_type( chemical::LOWER_TERMINUS_VARIANT ) );
		TS_ASSERT( l_pose.residue_type( 1 ).is_l_rna() );
	}

	// syn/anti chi restraints: an L residue gets the mirrored table, so mirror images score the same.
	void test_chi_restraints_are_mirrored() {
		pose::Pose d_pose, l_pose;
		make_mirror_pair( d_pose, l_pose );
		utility::vector1< Real > energies( 36, 20.0 );
		for ( Size q = 24; q <= 36; ++q ) energies[ q ] = 0.0; // the syn table of add_syn_chi_constraints
		scoring::ScoreFunction sfxn;
		sfxn.set_weight( scoring::rna_torsion, 1.0 );
		sfxn.set_weight( scoring::atom_pair_constraint, 1.0 ); // activates constraint scoring (the chi restraints are typed rna_torsion)
		Real const d_before = sfxn( d_pose ), l_before = sfxn( l_pose );
		for ( pose::Pose * pose : { &d_pose, &l_pose } ) {
			pose::full_model_info::make_sure_full_model_info_is_setup( *pose );
			pose::rna::add_chi_constraints_d_and_l( *pose, energies, utility::tools::make_vector1< Size >( 3, 7 ) );
		}
		Real const d_after = sfxn( d_pose ), l_after = sfxn( l_pose );
		TS_ASSERT( std::abs( d_after - d_before ) > 1.0 ); // the restraints do something
		TS_ASSERT_DELTA( l_after - l_before, d_after - d_before, 1.0e-3 );
	}

};
