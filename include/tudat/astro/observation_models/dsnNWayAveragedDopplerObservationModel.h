/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          T. Moyer (2000), Formulation for Observed and Computed Values of Deep Space Network Data Types for Navigation,
 *              DEEP SPACE COMMUNICATIONS AND NAVIGATION SERIES, JPL/NASA
 */

#ifndef TUDAT_DSNNWAYAVERAGEDDOPPLEROBSERVATIONMODEL_H
#define TUDAT_DSNNWAYAVERAGEDDOPPLEROBSERVATIONMODEL_H

#include <cmath>
#include <iostream>
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/simulation/simulation.h"

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/observationFrequencies.h"
#include "tudat/astro/observation_models/nWayRangeObservationModel.h"
#include "tudat/astro/observation_models/transmissionFrequencyInterface.h"
#include "tudat/interface/sofa/sofaTimeConversions.h"

namespace tudat
{

namespace observation_models
{

/*! Calculate the scaling factor for computing partials via DifferencedObservablePartial.
 *
 * Calculate the scaling factor for computing partials via DifferencedObservablePartial, for DSN n-way Doppler
 * observations. The scaling factor are selected according to eq. 13-59 of Moyer(2000).
 * The scaling factor depends on whether it is the first or the second partial being calculated. In
 * DifferencedObservablePartial, the first partial is multiplied by -1, hence here corresponds to the start of the
 * integration interval. The second partial is multiplied by +1, hence here corresponds to the end of the
 * integration interval.
 *
 * @param bodies System of bodies
 * @param linkEnds Map of the linkEnds defining the observation model
 * @param referenceLinkEnd Link end at which given time is valid, i.e. link end for which associated time
 *      is kept constant (to input value)
 * @param linkEndStates List of states at each link end during observation.
 * @param linkEndTimes List of times at each link end during observation.
 * @param ancillarySettings Observation ancillary simulation settings.
 * @param isFirstPartial Boolean indicating whether the scaling factor should be computed for the first (true) or
 *      second (false) partial.
 * @return Scaling factor
 */
inline double getDsnNWayAveragedDopplerScalingFactor(
        const std::function< double( std::vector< FrequencyBands > frequencyBands, double time ) > receivedFrequencyFunction,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< Eigen::Vector6d >& linkEndStates,
        const std::vector< double >& linkEndTimes,
        const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings,
        const bool isFirstPartial )
{
    double integrationTime;
    std::vector< FrequencyBands > frequencyBands;
    try
    {
        integrationTime = ancillarySettings->getAncilliaryDoubleData( doppler_integration_time );
        frequencyBands = convertDoubleVectorToFrequencyBands( ancillarySettings->getAncilliaryDoubleVectorData( frequency_bands ) );
    }
    catch( std::runtime_error& caughtException )
    {
        throw std::runtime_error( "Error when retrieving integration ancillary settings for DSN N-way averaged Doppler observable: " +
                                  std::string( caughtException.what( ) ) );
    }

    double transmissionTime;
    if( referenceLinkEnd == receiver )
    {
        if( isFirstPartial )
        {
            transmissionTime = linkEndTimes.at( 0 );
        }
        else
        {
            transmissionTime = linkEndTimes.at( 4 );
        }
    }
    //    else if ( referenceLinkEnd == transmitter )
    //    {
    //        if ( isFirstPartial )
    //        {
    //            transmissionTime = linkEndTimes.at( 3 );
    //        }
    //        else
    //        {
    //            transmissionTime = linkEndTimes.at( 7 );
    //        }
    //    }
    else
    {
        throw std::runtime_error( "Error when getting DSN N-way Doppler partials scaling factor: the selected reference link end (" +
                                  getLinkEndTypeString( referenceLinkEnd ) + ") is not valid." );
    }

    double frequency = receivedFrequencyFunction( frequencyBands, transmissionTime );

    // Moyer (2000), eq. 13-59
    return frequency / integrationTime / physical_constants::getSpeedOfLight< double >( );
}

/*! Compute the nodes and weights of an n-point Gauss-Legendre quadrature on [-1, 1].
 *
 * Nodes are the roots of the Legendre polynomial P_n, found by Newton iteration (in long double); the weights are
 * 2 / ( ( 1 - x^2 ) P_n'(x)^2 ).
 *
 * @param numberOfNodes Number of quadrature nodes (n >= 1)
 * @param nodes Quadrature nodes (output)
 * @param weights Quadrature weights (output)
 */
inline void computeGaussLegendreNodesAndWeights( const unsigned int numberOfNodes,
                                                 std::vector< long double >& nodes,
                                                 std::vector< long double >& weights )
{
    const long double pi = 3.141592653589793238462643383279502884L;
    nodes.assign( numberOfNodes, 0.0L );
    weights.assign( numberOfNodes, 0.0L );
    for( unsigned int i = 0; i < ( numberOfNodes + 1 ) / 2; i++ )
    {
        long double z = std::cos( pi * ( static_cast< long double >( i ) + 0.75L ) / ( static_cast< long double >( numberOfNodes ) + 0.5L ) );
        long double derivative = 1.0L;
        for( unsigned int iteration = 0; iteration < 100; iteration++ )
        {
            long double currentValue = 1.0L, previousValue = 0.0L;
            for( unsigned int j = 1; j <= numberOfNodes; j++ )
            {
                long double olderValue = previousValue;
                previousValue = currentValue;
                currentValue = ( ( 2.0L * j - 1.0L ) * z * previousValue - ( j - 1.0L ) * olderValue ) / j;
            }
            derivative = numberOfNodes * ( z * currentValue - previousValue ) / ( z * z - 1.0L );
            long double previousZ = z;
            z = previousZ - currentValue / derivative;
            if( std::fabs( z - previousZ ) < 1.0E-18L )
            {
                break;
            }
        }
        nodes.at( i ) = -z;
        nodes.at( numberOfNodes - 1 - i ) = z;
        weights.at( i ) = 2.0L / ( ( 1.0L - z * z ) * derivative * derivative );
        weights.at( numberOfNodes - 1 - i ) = weights.at( i );
    }
}

template< typename ObservationScalarType = double, typename TimeType = Time >
class DsnNWayAveragedDopplerObservationModel : public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;

    /*! Constructor.
     *
     * @param linkEnds Map of the linkEnds defining the observation model
     * @param arcStartObservationModel N-way range observation model associated with the start of the Doppler integration time.
     * @param arcEndObservationModel N-way range observation model associated with the end of the Doppler integration time.
     * @param bodyWithGroundStations Body object where the ground stations are located.
     * @param observationBiasCalculator Object for calculating (system-dependent) errors in the
     *  observable, i.e. deviations from the physically ideal observable between reference points (default none).
     */
    DsnNWayAveragedDopplerObservationModel(
            const LinkEnds& linkEnds,
            const std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcStartObservationModel,
            const std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcEndObservationModel,
            const std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittingFrequencyCalculator,
            const std::function< double( observation_models::FrequencyBands uplinkBand, observation_models::FrequencyBands downlinkBand ) >&
                    turnaroundRatio,
            const std::shared_ptr< ObservationBias< 1 > > observationBiasCalculator = nullptr,
            const std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > groundStationStates =
                    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > >( ),
            const bool subtractDopplerSignature = true ):
        ObservationModel< 1, ObservationScalarType, TimeType >( dsn_n_way_averaged_doppler, linkEnds, observationBiasCalculator ),
        arcStartObservationModel_( arcStartObservationModel ), arcEndObservationModel_( arcEndObservationModel ),
        numberOfLinkEnds_( linkEnds.size( ) ), transmittingFrequencyCalculator_( transmittingFrequencyCalculator ),
        turnaroundRatio_( turnaroundRatio ), stationStates_( groundStationStates ), subtractDopplerSignature_( subtractDopplerSignature )
    {
        if( !std::is_same< Time, TimeType >::value )
        {
            //            std::cerr<<
            //                    "Warning when defining DSN N-way averaged Doppler observation model: the selected time type "
            //                    "is not valid, using it would lead to large numerical errors."<<std::endl;
        }

        if( numberOfLinkEnds_ != 3 )
        {
            throw std::runtime_error(
                    "Error when defining DSN N-way averaged Doppler observation model: model allows exactly 3 link ends, " +
                    std::to_string( numberOfLinkEnds_ ) + "were selected." );
        }
        terrestrialTimeScaleConverter_ = earth_orientation::createDefaultTimeConverter( );

        computeGaussLegendreNodesAndWeights( numberOfQuadratureNodes_, quadratureNodes_, quadratureWeights_ );
    }

    /*! Set the functions used to compute the change of the signal path over the count without differencing light times.
     *
     * If set, the length of the transmission interval is computed from the displacement of each link end over the count
     * (see computeTransmissionUtcIntervalLength), which is free of the rounding errors of the (large) light times.
     * Otherwise, it is computed from the difference of the light times at the start and end of the count.
     *
     * @param linkEndBodyStateFunctions For each link end (in n-way order, transmitter first), function returning the
     *      state of the body on which the link end is located (in the global frame).
     * @param linkEndReferencePointOffsetFunctions For each link end, function returning the position of the link end
     *      w.r.t. the center of its body (in the global frame orientation), or nullptr if the link end is the body center.
     */
    void setLinkEndDisplacementFunctions(
            const std::vector< std::function< Eigen::Vector6d( const double ) > >& linkEndBodyStateFunctions,
            const std::vector< std::function< Eigen::Vector3d( const double ) > >& linkEndReferencePointOffsetFunctions )
    {
        if( linkEndBodyStateFunctions.size( ) != numberOfLinkEnds_ || linkEndReferencePointOffsetFunctions.size( ) != numberOfLinkEnds_ )
        {
            throw std::runtime_error( "Error when setting link end displacement functions of DSN N-way averaged Doppler model, size is "
                                      "inconsistent with number of link ends." );
        }
        linkEndBodyStateFunctions_ = linkEndBodyStateFunctions;
        linkEndReferencePointOffsetFunctions_ = linkEndReferencePointOffsetFunctions;
    }

    //! Destructor
    ~DsnNWayAveragedDopplerObservationModel( ) { }

    /*! Function to compute DSN n-way Doppler observation at given time.
     *
     * Function to compute DSN n-way Doppler observation at given time. Only implemented for receiver as the
     * linkEndAssociatedWithTime. Computes the observable according to section 13.3.2.2 of Moyer (2000).
     *
     * @param time Time at which observable is to be evaluated.
     * @param linkEndAssociatedWithTime Link end at which given time is valid, i.e. link end for which associated time
     *  is kept constant (to input value)
     * @param linkEndTimes List of times at each link end during observation.
     * @param linkEndStates List of states at each link end during observation.
     * @param ancillarySettings Observation ancillary simulation settings.
     * @return Observation value.
     */
    Eigen::Matrix< ObservationScalarType, 1, 1 > computeIdealObservationsWithLinkEndData(
            const TimeType time,
            const LinkEndType linkEndAssociatedWithTime,
            std::vector< double >& linkEndTimes,
            std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates,
            const std::shared_ptr< ObservationAncilliarySimulationSettings > ancillarySettings = nullptr )
    {
        // Check if selected reference link end is valid
        if( linkEndAssociatedWithTime != receiver )
        {
            throw std::runtime_error( "Error when computing DSN N-way Doppler observables: the selected reference link end (" +
                                      getLinkEndTypeString( linkEndAssociatedWithTime ) + ") is not valid." );
        }
        // Check if ancillary settings were provided
        if( ancillarySettings == nullptr )
        {
            throw std::runtime_error( "Error when simulating n-way DSN averaged Doppler observable; no ancillary settings found. " );
        }

        std::vector< double > arcStartLinkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > arcStartLinkEndStates;
        std::vector< double > arcEndLinkEndTimes;
        std::vector< Eigen::Matrix< double, 6, 1 > > arcEndLinkEndStates;

        TimeType integrationTime;
        ObservationScalarType referenceFrequency;
        std::vector< FrequencyBands > frequencyBands;
        FrequencyBands referenceUplinkBand;
        try
        {
            integrationTime = ancillarySettings->getAncilliaryDoubleData( doppler_integration_time );
            referenceFrequency = ancillarySettings->getAncilliaryDoubleData( doppler_reference_frequency );
            frequencyBands = convertDoubleVectorToFrequencyBands( ancillarySettings->getAncilliaryDoubleVectorData( frequency_bands ) );
            referenceUplinkBand =
                    convertDoubleToFrequencyBand( ancillarySettings->getAncilliaryDoubleData( reception_reference_frequency_band ) );
        }
        catch( std::runtime_error& caughtException )
        {
            throw std::runtime_error( "Error when retrieving ancillary settings for DSN N-way averaged Doppler observable: " +
                                      std::string( caughtException.what( ) ) );
        }

        if( frequencyBands.size( ) != numberOfLinkEnds_ - 1 )
        {
            throw std::runtime_error(
                    "Error when retrieving frequency bands ancillary settings for DSN N-way averaged Doppler observable: "
                    "size (" +
                    std::to_string( frequencyBands.size( ) ) + ") is inconsistent with number of links (" +
                    std::to_string( numberOfLinkEnds_ - 1 ) + ")." );
        }
        FrequencyBands uplinkBand = frequencyBands.at( 0 );
        FrequencyBands downlinkBand = frequencyBands.at( 1 );

        // Set approximate up- and down-link frequencies.
        ObservationScalarType currentTurnAroundRatio = static_cast< ObservationScalarType >( turnaroundRatio_( uplinkBand, downlinkBand ) );
        ObservationScalarType currentReferenceTurnAroundRatio =
                static_cast< ObservationScalarType >( turnaroundRatio_( referenceUplinkBand, downlinkBand ) );

        Eigen::Vector3d nominalReceivingStationState = ( stationStates_.count( receiver ) == 0 )
                ? Eigen::Vector3d::Zero( )
                : stationStates_.at( receiver )->getNominalCartesianPosition( );
        TimeType utcTime = terrestrialTimeScaleConverter_->getCurrentTime< TimeType >(
                basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, time, nominalReceivingStationState );

        TimeType receptionUtcStartTime = utcTime - integrationTime / 2.0;
        TimeType receptionUtcEndTime = utcTime + integrationTime / 2.0;

        TimeType receptionTdbStartTime = terrestrialTimeScaleConverter_->getCurrentTime< TimeType >(
                basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, receptionUtcStartTime, nominalReceivingStationState );
        TimeType receptionTdbEndTime = terrestrialTimeScaleConverter_->getCurrentTime< TimeType >(
                basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, receptionUtcEndTime, nominalReceivingStationState );

        Eigen::Vector3d nominalTransmittingStationState = ( stationStates_.count( transmitter ) == 0 )
                ? Eigen::Vector3d::Zero( )
                : stationStates_.at( transmitter )->getNominalCartesianPosition( );

        // Set frequencies for ionosphere/corona
        if( arcStartObservationModel_->getMultiLegLightTimeCalculator( )->doCorrectionsNeedFrequency( ) )
        {
            setTransmissionReceptionFrequencies( arcStartObservationModel_->getMultiLegLightTimeCalculator( ),
                                                 terrestrialTimeScaleConverter_,
                                                 transmittingFrequencyCalculator_,
                                                 receptionTdbStartTime,
                                                 linkEndAssociatedWithTime,
                                                 ancillarySettings,
                                                 currentTurnAroundRatio );
        }

        TimeType startLightTime =
                arcStartObservationModel_->computeIdealObservationsWithLinkEndData(
                        receptionTdbStartTime, linkEndAssociatedWithTime, arcStartLinkEndTimes, arcStartLinkEndStates, ancillarySettings )(
                        0, 0 ) /
                physical_constants::getSpeedOfLight< ObservationScalarType >( );

        // Light time and light-time correction of each leg at the start of the count
        std::vector< ObservationScalarType > arcStartLegLightTimes, arcStartLegCorrections;
        getLegLightTimesAndCorrections( arcStartObservationModel_, arcStartLegLightTimes, arcStartLegCorrections );

        // Set frequencies for ionosphere/corona
        if( arcEndObservationModel_->getMultiLegLightTimeCalculator( )->doCorrectionsNeedFrequency( ) )
        {
            setTransmissionReceptionFrequencies( arcEndObservationModel_->getMultiLegLightTimeCalculator( ),
                                                 terrestrialTimeScaleConverter_,
                                                 transmittingFrequencyCalculator_,
                                                 receptionTdbEndTime,
                                                 linkEndAssociatedWithTime,
                                                 ancillarySettings,
                                                 currentTurnAroundRatio );
        }
        TimeType endLightTime =
                arcEndObservationModel_->computeIdealObservationsWithLinkEndData(
                        receptionTdbEndTime, linkEndAssociatedWithTime, arcEndLinkEndTimes, arcEndLinkEndStates, ancillarySettings )( 0,
                                                                                                                                      0 ) /
                physical_constants::getSpeedOfLight< ObservationScalarType >( );

        // Light time and light-time correction of each leg at the end of the count
        std::vector< ObservationScalarType > arcEndLegLightTimes, arcEndLegCorrections;
        getLegLightTimesAndCorrections( arcEndObservationModel_, arcEndLegLightTimes, arcEndLegCorrections );

        // Moyer (2000), eqs. 13-52 and 13-53
        TimeType transmissionTdbStartTime = receptionTdbStartTime - startLightTime;
        TimeType transmissionTdbEndTime = receptionTdbEndTime - endLightTime;

        TimeType transmissionUtcStartTime = terrestrialTimeScaleConverter_->getCurrentTime< TimeType >(
                basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, transmissionTdbStartTime, nominalTransmittingStationState );
        TimeType transmissionUtcEndTime = terrestrialTimeScaleConverter_->getCurrentTime< TimeType >(
                basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, transmissionTdbEndTime, nominalTransmittingStationState );

        ObservationScalarType transmitterFrequencyIntegral;
        if( linkEndBodyStateFunctions_.size( ) != numberOfLinkEnds_ )
        {
            transmitterFrequencyIntegral =
                    transmittingFrequencyCalculator_->template getTemplatedFrequencyIntegral< ObservationScalarType, TimeType >(
                            transmissionUtcStartTime, transmissionUtcEndTime );
        }
        else
        {
            // The interval [transmissionUtcStartTime, transmissionUtcEndTime] above follows from differencing two light times
            // of up to several hours, whose rounding (~1E-12 s each, from the link end positions and the light time itself)
            // dominates the numerical error of this observable for short count times. Its length is instead computed from
            // quantities that are small by construction (see computeTransmissionUtcIntervalLength). The start time is kept:
            // an error in the placement of the interval only enters through the change of the frequency over that error.
            displacementCheckFailed_ = false;
            ObservationScalarType transmissionUtcIntervalLength = computeTransmissionUtcIntervalLength(
                    static_cast< ObservationScalarType >( integrationTime ),
                    receptionUtcStartTime,
                    receptionUtcEndTime,
                    receptionTdbStartTime,
                    receptionTdbEndTime,
                    transmissionTdbStartTime,
                    transmissionTdbEndTime,
                    transmissionUtcStartTime,
                    transmissionUtcEndTime,
                    nominalReceivingStationState,
                    nominalTransmittingStationState,
                    arcStartLinkEndTimes,
                    arcStartLinkEndStates,
                    arcStartLegLightTimes,
                    arcEndLegLightTimes,
                    arcStartLegCorrections,
                    arcEndLegCorrections );

            // If the velocity of a link-end body is not consistent with its position, use the legacy interval
            if( displacementCheckFailed_ )
            {
                if( !displacementCheckWarningGiven_ )
                {
                    std::cerr << "Warning in DSN N-way averaged Doppler model: ephemeris velocity of a link-end body is inconsistent "
                                 "with its position (displacement over the count differs by more than "
                              << maximumDisplacementDiscrepancy_
                              << " m); using the difference of the light times at the start and end of the count instead."
                              << std::endl;
                    displacementCheckWarningGiven_ = true;
                }
                transmissionUtcIntervalLength = static_cast< ObservationScalarType >( transmissionUtcEndTime - transmissionUtcStartTime );
            }

            // Integrate the frequency over an interval of this length. The end time is corrected to first order for its
            // rounding when represented as a TimeType (relevant if the seconds of the TimeType are a double).
            TimeType transmissionUtcIntervalEndTime = transmissionUtcStartTime + transmissionUtcIntervalLength;
            transmitterFrequencyIntegral =
                    transmittingFrequencyCalculator_->template getTemplatedFrequencyIntegral< ObservationScalarType, TimeType >(
                            transmissionUtcStartTime, transmissionUtcIntervalEndTime ) +
                    transmittingFrequencyCalculator_->template getTemplatedCurrentFrequency< ObservationScalarType, TimeType >(
                            transmissionUtcIntervalEndTime ) *
                            ( transmissionUtcIntervalLength -
                              static_cast< ObservationScalarType >( transmissionUtcIntervalEndTime - transmissionUtcStartTime ) );
        }

        // Moyer (2000), eq. 13-54
        Eigen::Matrix< ObservationScalarType, 1, 1 > observation =
                ( Eigen::Matrix< ObservationScalarType, 1, 1 >( ) << currentReferenceTurnAroundRatio * referenceFrequency +
                          ( subtractDopplerSignature_ ? mathematical_constants::getFloatingInteger< ObservationScalarType >( -1.0 )
                                                      : mathematical_constants::getFloatingInteger< ObservationScalarType >( 1.0 ) ) *
                                  currentTurnAroundRatio / static_cast< ObservationScalarType >( integrationTime ) *
                                  transmitterFrequencyIntegral )
                        .finished( );

        linkEndTimes.clear( );
        linkEndStates.clear( );
        linkEndTimes.resize( 4 * ( numberOfLinkEnds_ - 1 ) );
        linkEndStates.resize( 4 * ( numberOfLinkEnds_ - 1 ) );

        for( unsigned int i = 0; i < 2 * ( numberOfLinkEnds_ - 1 ); i++ )
        {
            linkEndTimes[ i ] = arcStartLinkEndTimes[ i ];
            linkEndTimes[ i + 2 * ( numberOfLinkEnds_ - 1 ) ] = arcEndLinkEndTimes[ i ];

            linkEndStates[ i ] = arcStartLinkEndStates[ i ];
            linkEndStates[ i + 2 * ( numberOfLinkEnds_ - 1 ) ] = arcEndLinkEndStates[ i ];
        }

        return observation;
    }

    // Function to retrieve the arc end observation model
    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > getArcEndObservationModel( )
    {
        return arcEndObservationModel_;
    }

    // Function to retrieve the arc start observation model
    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > getArcStartObservationModel( )
    {
        return arcStartObservationModel_;
    }

private:
    // Retrieve the light time and light-time correction of each leg of the last solution of the given n-way range model
    void getLegLightTimesAndCorrections( const std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > nWayRangeModel,
                                         std::vector< ObservationScalarType >& legLightTimes,
                                         std::vector< ObservationScalarType >& legCorrections )
    {
        legLightTimes.clear( );
        legCorrections.clear( );
        for( auto legCalculator: nWayRangeModel->getMultiLegLightTimeCalculator( )->getLightTimeCalculators( ) )
        {
            legLightTimes.push_back( legCalculator->getCurrentIdealLightTime( ) + legCalculator->getCurrentLightTimeCorrection( ) );
            legCorrections.push_back( legCalculator->getCurrentLightTimeCorrection( ) );
        }
    }

    // Displacement of the body (center) of the given link end over [startTime, startTime + intervalLength], by composite
    // Gauss-Legendre quadrature of its velocity
    Eigen::Matrix< ObservationScalarType, 3, 1 > computeBodyDisplacement( const unsigned int linkEndIndex,
                                                                        const double startTime,
                                                                        const ObservationScalarType intervalLength )
    {
        unsigned int numberOfSubintervals = 1;
        if( maximumQuadratureSubintervalLength_ > 0.0 )
        {
            numberOfSubintervals = std::max(
                    1, static_cast< int >( std::ceil( std::fabs( static_cast< double >( intervalLength ) ) / maximumQuadratureSubintervalLength_ ) ) );
        }
        ObservationScalarType subintervalLength = intervalLength / static_cast< ObservationScalarType >( numberOfSubintervals );

        Eigen::Matrix< ObservationScalarType, 3, 1 > displacement = Eigen::Matrix< ObservationScalarType, 3, 1 >::Zero( );
        for( unsigned int subinterval = 0; subinterval < numberOfSubintervals; subinterval++ )
        {
            ObservationScalarType subintervalCenter = ( static_cast< ObservationScalarType >( subinterval ) + 0.5 ) * subintervalLength;
            for( unsigned int node = 0; node < quadratureNodes_.size( ); node++ )
            {
                double nodeTime = startTime +
                        static_cast< double >( subintervalCenter +
                                               static_cast< ObservationScalarType >( quadratureNodes_.at( node ) ) * subintervalLength / 2.0 );
                displacement += static_cast< ObservationScalarType >( quadratureWeights_.at( node ) ) * subintervalLength / 2.0 *
                        linkEndBodyStateFunctions_.at( linkEndIndex )( nodeTime ).segment( 3, 3 ).template cast< ObservationScalarType >( );
            }
        }

        // Check against the difference of the positions (rounded, but free of any inconsistency between the velocity and
        // position of the ephemeris, e.g. a custom ephemeris with zero or approximate velocity)
        double endTime = startTime + static_cast< double >( intervalLength );
        Eigen::Vector6d startState = linkEndBodyStateFunctions_.at( linkEndIndex )( startTime );
        Eigen::Vector6d endState = linkEndBodyStateFunctions_.at( linkEndIndex )( endTime );
        Eigen::Vector3d positionDifference = endState.segment( 0, 3 ) - startState.segment( 0, 3 ) +
                endState.segment( 3, 3 ) * ( static_cast< double >( intervalLength ) - ( endTime - startTime ) );
        double displacementDiscrepancy = ( displacement.template cast< double >( ) - positionDifference ).norm( );
        if( displacementDiscrepancy > maximumDisplacementDiscrepancy_ )
        {
            displacementCheckFailed_ = true;
        }
        return displacement;
    }

    // Change of the position of the given link end w.r.t. its body over [startTime, startTime + intervalLength] (zero if the
    // link end is the body center). The positions are small, so they can be differenced directly; the end time is corrected
    // (to first order) for its rounding when represented as a double.
    Eigen::Matrix< ObservationScalarType, 3, 1 > computeReferencePointDisplacement( const unsigned int linkEndIndex,
                                                                                  const double startTime,
                                                                                  const ObservationScalarType intervalLength )
    {
        if( linkEndReferencePointOffsetFunctions_.at( linkEndIndex ) == nullptr )
        {
            return Eigen::Matrix< ObservationScalarType, 3, 1 >::Zero( );
        }
        double endTime = startTime + static_cast< double >( intervalLength );
        Eigen::Matrix< ObservationScalarType, 3, 1 > displacement =
                ( linkEndReferencePointOffsetFunctions_.at( linkEndIndex )( endTime ) -
                  linkEndReferencePointOffsetFunctions_.at( linkEndIndex )( startTime ) )
                        .template cast< ObservationScalarType >( );
        ObservationScalarType representedLength = static_cast< ObservationScalarType >( endTime - startTime );
        if( representedLength != 0.0 )
        {
            displacement += displacement / representedLength * ( intervalLength - representedLength );
        }
        return displacement;
    }

    /*! Compute the length (in UTC at the transmitting station) of the transmission interval of the count.
     *
     * The length is L_T = L_R - sum_j ( Delta rho_j / c + Delta C_j ) - Delta( TDB - UTC )_T, with
     *  - L_R = T + Delta( TDB - UTC )_R the length of the reception interval in TDB (T the count time in UTC),
     *  - Delta rho_j the change of the Euclidean length of leg j over the count, from the change of the leg vector,
     *    Delta r_j = Delta x_R - Delta x_T. The displacement of each link end over its own interval is the displacement of
     *    its body (quadrature of the body velocity) plus the change of its position w.r.t. the body (a small vector);
     *    Delta rho_j = ( 2 r_j . Delta r_j + |Delta r_j|^2 ) / ( |r_j + Delta r_j| + |r_j| ), with r_j the leg vector at
     *    the start of the count,
     *  - Delta C_j the change of the light-time correction of leg j between the solutions at the start and end of the
     *    count (retransmission delays are constant and do not contribute),
     *  - Delta( TDB - UTC ) the change of TDB - TT over the interval at the station (plus any leap second).
     * None of these terms requires differencing the two (large) light times at the start and end of the count; the
     * interval of each link end uses the (rounded) light times, whose error only enters multiplied by v/c.
     */
    ObservationScalarType computeTransmissionUtcIntervalLength( const ObservationScalarType integrationTime,
                                                                const TimeType receptionUtcStartTime,
                                                                const TimeType receptionUtcEndTime,
                                                                const TimeType receptionTdbStartTime,
                                                                const TimeType receptionTdbEndTime,
                                                                const TimeType transmissionTdbStartTime,
                                                                const TimeType transmissionTdbEndTime,
                                                                const TimeType transmissionUtcStartTime,
                                                                const TimeType transmissionUtcEndTime,
                                                                const Eigen::Vector3d& nominalReceivingStationState,
                                                                const Eigen::Vector3d& nominalTransmittingStationState,
                                                                const std::vector< double >& arcStartLinkEndTimes,
                                                                const std::vector< Eigen::Matrix< double, 6, 1 > >& arcStartLinkEndStates,
                                                                const std::vector< ObservationScalarType >& arcStartLegLightTimes,
                                                                const std::vector< ObservationScalarType >& arcEndLegLightTimes,
                                                                const std::vector< ObservationScalarType >& arcStartLegCorrections,
                                                                const std::vector< ObservationScalarType >& arcEndLegCorrections )
    {
        // Length of the reception interval in TDB: change of TDB - TT over the interval (plus the change in leap seconds,
        // an integer obtained from the direct difference of the converted times)
        double receptionTdbMinusTtChange =
                sofa_interface::getTDBminusTT( static_cast< double >( receptionTdbEndTime ), nominalReceivingStationState ) -
                sofa_interface::getTDBminusTT( static_cast< double >( receptionTdbStartTime ), nominalReceivingStationState );
        double receptionLeapSecondChange = std::round(
                static_cast< double >( ( receptionTdbEndTime - receptionTdbStartTime ) - ( receptionUtcEndTime - receptionUtcStartTime ) ) -
                receptionTdbMinusTtChange );
        ObservationScalarType receptionTdbIntervalLength = integrationTime + static_cast< ObservationScalarType >( receptionTdbMinusTtChange ) +
                static_cast< ObservationScalarType >( receptionLeapSecondChange );

        // Change of each leg's light time, moving from the receiver to the transmitter
        const unsigned int numberOfLegs = arcStartLegLightTimes.size( );
        ObservationScalarType legIntervalLength = receptionTdbIntervalLength;
        ObservationScalarType totalLightTimeChange = mathematical_constants::getFloatingInteger< ObservationScalarType >( 0 );
        for( int leg = static_cast< int >( numberOfLegs ) - 1; leg >= 0; leg-- )
        {
            // Link end intervals of this leg: the reception interval is that of the transmission of the next leg (or of the
            // count), the transmission interval uses the change of the (rounded) light time of the leg
            ObservationScalarType legReceptionIntervalLength = legIntervalLength;
            ObservationScalarType legTransmissionIntervalLength =
                    legReceptionIntervalLength - ( arcEndLegLightTimes.at( leg ) - arcStartLegLightTimes.at( leg ) );

            // Change of the leg vector over the count
            Eigen::Matrix< ObservationScalarType, 3, 1 > legVectorChange =
                    ( computeBodyDisplacement( leg + 1, arcStartLinkEndTimes.at( 2 * leg + 1 ), legReceptionIntervalLength ) +
                      computeReferencePointDisplacement( leg + 1, arcStartLinkEndTimes.at( 2 * leg + 1 ), legReceptionIntervalLength ) ) -
                    ( computeBodyDisplacement( leg, arcStartLinkEndTimes.at( 2 * leg ), legTransmissionIntervalLength ) +
                      computeReferencePointDisplacement( leg, arcStartLinkEndTimes.at( 2 * leg ), legTransmissionIntervalLength ) );

            // Change of the Euclidean leg length, computed without differencing the (large) lengths themselves
            Eigen::Matrix< ObservationScalarType, 3, 1 > legVector =
                    ( arcStartLinkEndStates.at( 2 * leg + 1 ).segment( 0, 3 ) - arcStartLinkEndStates.at( 2 * leg ).segment( 0, 3 ) )
                            .template cast< ObservationScalarType >( );
            ObservationScalarType legLengthChange =
                    ( 2.0 * legVector.dot( legVectorChange ) + legVectorChange.dot( legVectorChange ) ) /
                    ( ( legVector + legVectorChange ).norm( ) + legVector.norm( ) );

            ObservationScalarType legLightTimeChange = legLengthChange / physical_constants::getSpeedOfLight< ObservationScalarType >( ) +
                    ( arcEndLegCorrections.at( leg ) - arcStartLegCorrections.at( leg ) );
            totalLightTimeChange += legLightTimeChange;

            // Transmission interval of this leg is the reception interval of the previous one (retransmission delays constant)
            legIntervalLength = legReceptionIntervalLength - legLightTimeChange;
        }

        // Length of the transmission interval in TDB
        ObservationScalarType transmissionTdbIntervalLength = receptionTdbIntervalLength - totalLightTimeChange;

        // Length of the transmission interval in UTC
        double transmissionTdbMinusTtChange =
                sofa_interface::getTDBminusTT(
                        static_cast< double >( transmissionTdbStartTime ) + static_cast< double >( transmissionTdbIntervalLength ),
                        nominalTransmittingStationState ) -
                sofa_interface::getTDBminusTT( static_cast< double >( transmissionTdbStartTime ), nominalTransmittingStationState );
        double transmissionLeapSecondChange =
                std::round( static_cast< double >( ( transmissionTdbEndTime - transmissionTdbStartTime ) -
                                                   ( transmissionUtcEndTime - transmissionUtcStartTime ) ) -
                            transmissionTdbMinusTtChange );
        return transmissionTdbIntervalLength - static_cast< ObservationScalarType >( transmissionTdbMinusTtChange ) -
                static_cast< ObservationScalarType >( transmissionLeapSecondChange );
    }

    // For each link end (n-way order), function returning the state of its body in the global frame
    std::vector< std::function< Eigen::Vector6d( const double ) > > linkEndBodyStateFunctions_;

    // For each link end (n-way order), function returning its position w.r.t. its body (nullptr for a body center)
    std::vector< std::function< Eigen::Vector3d( const double ) > > linkEndReferencePointOffsetFunctions_;

    // Number of Gauss-Legendre nodes (per subinterval) used to compute the displacement of the link-end bodies
    unsigned int numberOfQuadratureNodes_ = 8;

    // Maximum length of a quadrature subinterval (s); longer intervals use a composite rule
    double maximumQuadratureSubintervalLength_ = 60.0;

    // Maximum difference (m) between the displacement of a link-end body from its velocity and from its positions
    double maximumDisplacementDiscrepancy_ = 1.0E-2;

    // Whether the displacement check failed for the current observation, and whether a warning was given
    bool displacementCheckFailed_ = false;
    bool displacementCheckWarningGiven_ = false;

    // Gauss-Legendre nodes and weights on [-1, 1]
    std::vector< long double > quadratureNodes_;
    std::vector< long double > quadratureWeights_;

    // N-way range observation model associated with the start of the Doppler integration time.
    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcStartObservationModel_;

    // N-way range observation model associated with the end of the Doppler integration time.
    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcEndObservationModel_;

    // Number of link ends
    unsigned int numberOfLinkEnds_;

    // Object returning the transmitted frequency as the transmitting link end
    std::shared_ptr< ground_stations::StationFrequencyInterpolator > transmittingFrequencyCalculator_;

    // Function returning the turnaround ratio for given uplink and downlink bands
    std::function< double( FrequencyBands uplinkBand, FrequencyBands downlinkBand ) > turnaroundRatio_;

    std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter > terrestrialTimeScaleConverter_;

    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > stationStates_;

    bool subtractDopplerSignature_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_DSNNWAYAVERAGEDDOPPLEROBSERVATIONMODEL_H
