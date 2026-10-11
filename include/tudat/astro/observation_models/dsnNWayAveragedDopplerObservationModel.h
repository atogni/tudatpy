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
#include <functional>
#include <stdexcept>
#include <string>
#include <vector>

#include "tudat/astro/ground_stations/groundStationState.h"
#include "tudat/astro/ground_stations/transmittingFrequencies.h"
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
template< typename ObservationScalarType >
inline double getDsnNWayAveragedDopplerScalingFactor(
        const std::function< double( std::vector< FrequencyBands > frequencyBands, double time ) > receivedFrequencyFunction,
        const bool subtractDopplerSignature,
        const observation_models::LinkEndType referenceLinkEnd,
        const std::vector< Eigen::Vector6d >& linkEndStates,
        const std::vector< double >& linkEndTimes,
        const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings,
        const bool isFirstPartial )
{
    double integrationTime;
    std::vector< FrequencyBands > frequencyBands;
    try
    {
        integrationTime = ancillarySettings->getAncillaryDoubleData( doppler_integration_time );
        frequencyBands = convertDoubleVectorToFrequencyBands( ancillarySettings->getAncillaryDoubleVectorData( frequency_bands ) );
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
    return ( subtractDopplerSignature ? mathematical_constants::getFloatingInteger< ObservationScalarType >( 1.0 )
                                      : mathematical_constants::getFloatingInteger< ObservationScalarType >( -1.0 ) ) *
            frequency / integrationTime / physical_constants::getSpeedOfLight< double >( );
}

template< typename ObservationScalarType = double, typename TimeType = Time >
class DsnNWayAveragedDopplerObservationModel : public ObservationModel< 1, ObservationScalarType, TimeType >
{
public:
    typedef Eigen::Matrix< ObservationScalarType, 6, 1 > StateType;

    using ObservationModel< 1, ObservationScalarType, TimeType >::turnaroundRatio_;
    using ObservationModel< 1, ObservationScalarType, TimeType >::timeScaleConverter_;
    using ObservationModel< 1, ObservationScalarType, TimeType >::frequencyInterpolator_;

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
        ObservationModel< 1, ObservationScalarType, TimeType >(
                dsn_n_way_averaged_doppler,
                linkEnds,
                observationBiasCalculator,
                std::vector< std::shared_ptr< FullLinkLightTimeCalculator< ObservationScalarType, TimeType > > >{
                        arcStartObservationModel->getFullLinkLightTimeCalculator( ),
                        arcEndObservationModel->getFullLinkLightTimeCalculator( ) } ),
        arcStartObservationModel_( arcStartObservationModel ), arcEndObservationModel_( arcEndObservationModel ),
        numberOfLinkEnds_( linkEnds.size( ) ), stationStates_( groundStationStates ), subtractDopplerSignature_( subtractDopplerSignature )
    {
        this->setFrequencyInterpolatorAndTurnaroundRatio( transmittingFrequencyCalculator, turnaroundRatio );
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
    }

    //! Destructor
    ~DsnNWayAveragedDopplerObservationModel( ) {}

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
            const std::shared_ptr< ObservationAncillarySimulationSettings > ancillarySettings = nullptr ) override
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
            integrationTime = ancillarySettings->getAncillaryDoubleData( doppler_integration_time );
            referenceFrequency = ancillarySettings->getAncillaryDoubleData( doppler_reference_frequency );
            frequencyBands = convertDoubleVectorToFrequencyBands( ancillarySettings->getAncillaryDoubleVectorData( frequency_bands ) );
            referenceUplinkBand =
                    convertDoubleToFrequencyBand( ancillarySettings->getAncillaryDoubleData( reception_reference_frequency_band ) );
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
        TimeType utcTime = timeScaleConverter_->template getCurrentTime< TimeType >(
                basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, time, nominalReceivingStationState );

        TimeType receptionUtcStartTime = utcTime - integrationTime / 2.0;
        TimeType receptionUtcEndTime = utcTime + integrationTime / 2.0;

        TimeType receptionTdbStartTime = timeScaleConverter_->template getCurrentTime< TimeType >(
                basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, receptionUtcStartTime, nominalReceivingStationState );
        TimeType receptionTdbEndTime = timeScaleConverter_->template getCurrentTime< TimeType >(
                basic_astrodynamics::utc_scale, basic_astrodynamics::tdb_scale, receptionUtcEndTime, nominalReceivingStationState );

        Eigen::Vector3d nominalTransmittingStationState = ( stationStates_.count( transmitter ) == 0 )
                ? Eigen::Vector3d::Zero( )
                : stationStates_.at( transmitter )->getNominalCartesianPosition( );

        // Set frequencies for ionosphere/corona
        if( arcStartObservationModel_->getFullLinkLightTimeCalculator( )->doCorrectionsNeedFrequency( ) )
        {
            setTransmissionReceptionFrequencies( arcStartObservationModel_->getFullLinkLightTimeCalculator( ),
                                                 timeScaleConverter_,
                                                 frequencyInterpolator_,
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

        // Set frequencies for ionosphere/corona
        if( arcEndObservationModel_->getFullLinkLightTimeCalculator( )->doCorrectionsNeedFrequency( ) )
        {
            setTransmissionReceptionFrequencies( arcEndObservationModel_->getFullLinkLightTimeCalculator( ),
                                                 timeScaleConverter_,
                                                 frequencyInterpolator_,
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

        // Moyer (2000), eqs. 13-52 and 13-53
        TimeType transmissionTdbStartTime = receptionTdbStartTime - startLightTime;
        TimeType transmissionTdbEndTime = receptionTdbEndTime - endLightTime;

        TimeType transmissionUtcStartTime = timeScaleConverter_->template getCurrentTime< TimeType >(
                basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, transmissionTdbStartTime, nominalTransmittingStationState );
        TimeType transmissionUtcEndTime = timeScaleConverter_->template getCurrentTime< TimeType >(
                basic_astrodynamics::tdb_scale, basic_astrodynamics::utc_scale, transmissionTdbEndTime, nominalTransmittingStationState );

        ObservationScalarType transmitterFrequencyIntegral;
        ObservationScalarType transmissionUtcIntervalLength = computeTransmissionUtcIntervalLength( integrationTime,
                                                                                                    receptionTdbStartTime,
                                                                                                    receptionTdbEndTime,
                                                                                                    transmissionTdbStartTime,
                                                                                                    transmissionTdbEndTime,
                                                                                                    transmissionUtcStartTime,
                                                                                                    transmissionUtcEndTime,
                                                                                                    nominalReceivingStationState,
                                                                                                    nominalTransmittingStationState,
                                                                                                    arcStartLinkEndTimes,
                                                                                                    arcStartLinkEndStates );
        if( std::isnan( transmissionUtcIntervalLength ) )
        {
            transmitterFrequencyIntegral =
                    frequencyInterpolator_->template getTemplatedFrequencyIntegral< ObservationScalarType, TimeType >(
                            transmissionUtcStartTime, transmissionUtcEndTime );
        }
        else
        {
            // Integrate over an interval of the computed length; the end time is corrected to first order for its rounding
            // when represented as a TimeType
            TimeType transmissionUtcIntervalEndTime = transmissionUtcStartTime + transmissionUtcIntervalLength;
            transmitterFrequencyIntegral =
                    frequencyInterpolator_->template getTemplatedFrequencyIntegral< ObservationScalarType, TimeType >(
                            transmissionUtcStartTime, transmissionUtcIntervalEndTime ) +
                    frequencyInterpolator_->template getTemplatedCurrentFrequency< ObservationScalarType, TimeType >(
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

    bool getSubtractDopplerSignature( )
    {
        return subtractDopplerSignature_;
    }

    /*! Set the link-end motion (per link end, transmitter first: body state in the global frame; position w.r.t. the body, or
     * nullptr for a body center), used to compute the light-time change over the count from link-end displacements instead
     * of differencing the two light times, whose rounding dominates the numerical error at short count times.
     */
    void setLinkEndMotionFunctions( const std::vector< std::function< Eigen::Vector6d( const double ) > >& linkEndBodyStateFunctions,
                                    const std::vector< std::function< Eigen::Vector3d( const double ) > >& linkEndStationPositionFunctions )
    {
        linkEndBodyStateFunctions_ = linkEndBodyStateFunctions;
        linkEndStationPositionFunctions_ = linkEndStationPositionFunctions;
    }

    std::map< std::pair< LinkEndType, LinkEndType >, std::vector< std::shared_ptr< LightTimeCalculatorBase > > >
    getLegLightTimeCalculators( ) const override
    {
        std::map< std::pair< LinkEndType, LinkEndType >, std::vector< std::shared_ptr< LightTimeCalculatorBase > > > legMap =
                arcStartObservationModel_->getLegLightTimeCalculators( );
        const std::map< std::pair< LinkEndType, LinkEndType >, std::vector< std::shared_ptr< LightTimeCalculatorBase > > > endLegMap =
                arcEndObservationModel_->getLegLightTimeCalculators( );
        for( const auto& endLegEntry : endLegMap )
        {
            legMap[ endLegEntry.first ].insert(
                    legMap[ endLegEntry.first ].end( ), endLegEntry.second.begin( ), endLegEntry.second.end( ) );
        }
        return legMap;
    }

private:
    /*! Length (UTC at the transmitter) of the transmission interval, Moyer (2000) eq. 13-51, without differencing the two light
     * times (NaN if no link-end motion is set or a body's velocity is inconsistent with its positions: the caller then
     * differences them). Light-time change of leg j: Delta rho_j / c + Delta C_j, with
     * Delta rho_j = ( 2 r_j . Delta r_j + |Delta r_j|^2 ) / ( |r_j + Delta r_j| + |r_j| ) and Delta r_j from the link-end
     * displacements (body velocity integrated by Gauss-Legendre quadrature, plus the change of the position w.r.t. the body).
     */
    ObservationScalarType computeTransmissionUtcIntervalLength( const TimeType integrationTime,
                                                                const TimeType receptionTdbStartTime,
                                                                const TimeType receptionTdbEndTime,
                                                                const TimeType transmissionTdbStartTime,
                                                                const TimeType transmissionTdbEndTime,
                                                                const TimeType transmissionUtcStartTime,
                                                                const TimeType transmissionUtcEndTime,
                                                                const Eigen::Vector3d& receivingStationPosition,
                                                                const Eigen::Vector3d& transmittingStationPosition,
                                                                const std::vector< double >& arcStartLinkEndTimes,
                                                                const std::vector< Eigen::Matrix< double, 6, 1 > >& arcStartLinkEndStates )
    {
        if( linkEndBodyStateFunctions_.size( ) != 3 )
        {
            return TUDAT_NAN;
        }
        // 5-point Gauss-Legendre quadrature on [-1, 1]
        const double nodes[ 5 ] = { -0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640 };
        const double weights[ 5 ] = { 0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665, 0.2369268850561891 };

        // Reception interval in TDB: count time plus the change of TDB - TT (and of leap seconds) at the receiver
        double receptionTdbMinusTtChange =
                sofa_interface::getTDBminusTT( static_cast< double >( receptionTdbEndTime ), receivingStationPosition ) -
                sofa_interface::getTDBminusTT( static_cast< double >( receptionTdbStartTime ), receivingStationPosition );
        double receptionLeapSecondChange = std::round(
                static_cast< double >( ( receptionTdbEndTime - receptionTdbStartTime ) - integrationTime ) - receptionTdbMinusTtChange );
        ObservationScalarType receptionIntervalLength = static_cast< ObservationScalarType >( integrationTime ) +
                static_cast< ObservationScalarType >( receptionTdbMinusTtChange + receptionLeapSecondChange );

        // Interval of each link end (transmitter, retransmitter, receiver)
        std::vector< std::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType > > > startLegs =
                arcStartObservationModel_->getLightTimeCalculators( );
        std::vector< std::shared_ptr< LightTimeCalculator< ObservationScalarType, TimeType > > > endLegs =
                arcEndObservationModel_->getLightTimeCalculators( );
        std::vector< ObservationScalarType > intervalLengths = { 0.0, 0.0, receptionIntervalLength };
        for( int leg = 1; leg >= 0; leg-- )
        {
            intervalLengths.at( leg ) = intervalLengths.at( leg + 1 ) -
                    ( endLegs.at( leg )->getCurrentIdealLightTime( ) + endLegs.at( leg )->getCurrentLightTimeCorrection( ) -
                      startLegs.at( leg )->getCurrentIdealLightTime( ) - startLegs.at( leg )->getCurrentLightTimeCorrection( ) );
        }

        // Displacement of each link end over its interval
        std::vector< Eigen::Matrix< ObservationScalarType, 3, 1 > > displacements( 3 );
        for( unsigned int i = 0; i < 3; i++ )
        {
            double startTime = arcStartLinkEndTimes.at( i == 0 ? 0 : 2 * i - 1 );
            double endTime = startTime + static_cast< double >( intervalLengths.at( i ) );
            int numberOfSubintervals = std::max( 1, static_cast< int >( std::ceil( std::fabs( endTime - startTime ) / 60.0 ) ) );
            ObservationScalarType subintervalLength =
                    intervalLengths.at( i ) / static_cast< ObservationScalarType >( numberOfSubintervals );
            displacements.at( i ).setZero( );
            for( int j = 0; j < numberOfSubintervals; j++ )
            {
                for( unsigned int k = 0; k < 5; k++ )
                {
                    double nodeTime = startTime + static_cast< double >( ( j + 0.5 + 0.5 * nodes[ k ] ) * subintervalLength );
                    displacements.at( i ) += static_cast< ObservationScalarType >( 0.5 * weights[ k ] ) * subintervalLength *
                            linkEndBodyStateFunctions_.at( i )( nodeTime ).segment( 3, 3 ).template cast< ObservationScalarType >( );
                }
            }

            // Guard against an ephemeris whose velocity is inconsistent with its positions
            Eigen::Vector3d positionChange = linkEndBodyStateFunctions_.at( i )( endTime ).segment( 0, 3 ) -
                    linkEndBodyStateFunctions_.at( i )( startTime ).segment( 0, 3 );
            if( ( displacements.at( i ).template cast< double >( ) - positionChange ).norm( ) > 1.0E-2 )
            {
                return TUDAT_NAN;
            }

            // Change of the position w.r.t. the body (small, so differenced directly; scaled for the rounding of endTime)
            if( linkEndStationPositionFunctions_.at( i ) != nullptr && endTime != startTime )
            {
                displacements.at( i ) +=
                        ( ( linkEndStationPositionFunctions_.at( i )( endTime ) - linkEndStationPositionFunctions_.at( i )( startTime ) ) *
                          ( static_cast< double >( intervalLengths.at( i ) ) / ( endTime - startTime ) ) )
                                .template cast< ObservationScalarType >( );
            }
        }

        // Change of the total light time
        ObservationScalarType lightTimeChange = 0.0;
        for( unsigned int leg = 0; leg < 2; leg++ )
        {
            Eigen::Matrix< ObservationScalarType, 3, 1 > legVector =
                    ( arcStartLinkEndStates.at( 2 * leg + 1 ).segment( 0, 3 ) - arcStartLinkEndStates.at( 2 * leg ).segment( 0, 3 ) )
                            .template cast< ObservationScalarType >( );
            Eigen::Matrix< ObservationScalarType, 3, 1 > legVectorChange = displacements.at( leg + 1 ) - displacements.at( leg );
            lightTimeChange += ( 2.0 * legVector.dot( legVectorChange ) + legVectorChange.squaredNorm( ) ) /
                            ( ( legVector + legVectorChange ).norm( ) + legVector.norm( ) ) /
                            physical_constants::getSpeedOfLight< ObservationScalarType >( ) +
                    ( endLegs.at( leg )->getCurrentLightTimeCorrection( ) - startLegs.at( leg )->getCurrentLightTimeCorrection( ) );
        }

        // Transmission interval in TDB, then in UTC at the transmitter
        ObservationScalarType transmissionTdbIntervalLength = receptionIntervalLength - lightTimeChange;
        double transmissionTdbMinusTtChange = sofa_interface::getTDBminusTT( static_cast< double >( transmissionTdbStartTime ) +
                                                                                     static_cast< double >( transmissionTdbIntervalLength ),
                                                                             transmittingStationPosition ) -
                sofa_interface::getTDBminusTT( static_cast< double >( transmissionTdbStartTime ), transmittingStationPosition );
        double transmissionLeapSecondChange = std::round( static_cast< double >( ( transmissionTdbEndTime - transmissionTdbStartTime ) -
                                                                                 ( transmissionUtcEndTime - transmissionUtcStartTime ) ) -
                                                          transmissionTdbMinusTtChange );
        return transmissionTdbIntervalLength -
                static_cast< ObservationScalarType >( transmissionTdbMinusTtChange + transmissionLeapSecondChange );
    }

    // Link-end motion (see setLinkEndMotionFunctions)
    std::vector< std::function< Eigen::Vector6d( const double ) > > linkEndBodyStateFunctions_;
    std::vector< std::function< Eigen::Vector3d( const double ) > > linkEndStationPositionFunctions_;

    // N-way range observation model associated with the start of the Doppler integration time.
    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcStartObservationModel_;

    // N-way range observation model associated with the end of the Doppler integration time.
    std::shared_ptr< NWayRangeObservationModel< ObservationScalarType, TimeType > > arcEndObservationModel_;

    // Number of link ends
    unsigned int numberOfLinkEnds_;

    std::map< LinkEndType, std::shared_ptr< ground_stations::GroundStationState > > stationStates_;

    bool subtractDopplerSignature_;
};

}  // namespace observation_models

}  // namespace tudat

#endif  // TUDAT_DSNNWAYAVERAGEDDOPPLEROBSERVATIONMODEL_H
