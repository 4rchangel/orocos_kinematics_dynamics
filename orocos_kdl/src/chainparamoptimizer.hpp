// Copyright  (C)  2022  Gabriel Pfeilschifter <gabriel dot pfeilschifter at posteo dot de>

// Version: 1.0
// Author: Gabriel Pfeilschifter <gabriel dot pfeilschifter at posteo dot de>
// Maintainer: Gabriel Pfeilschifter <gabriel dot pfeilschifter at posteo dot de>
// URL: http://www.orocos.org/kdl

// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.

// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.

// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

#ifndef KDL_CHAINPARAMOPTIMIZER_HPP
#define KDL_CHAINPARAMOPTIMIZER_HPP

#include <memory>
#include "frames.hpp"
#include "jntarray.hpp"
#include "chain.hpp"
#include "solveri.hpp"

namespace KDL
{
    class SegmentDH : public Segment
	{
	public:
		SegmentDH(const std::string& name, const Joint& joint, double a,double alpha,double d,double theta, const RigidBodyInertia& _I = RigidBodyInertia::Zero());
		SegmentDH(const SegmentDH& other);
		SegmentDH& operator=(const SegmentDH& arg);

		double GetDHParamA() const {return dhParamA;}
		double GetDHParamD() const {return dhParamD;}
		double GetDHParamAlpha() const {return dhParamAlpha;}
		double GetDHParamTeta() const {return dhParamTheta;}

		virtual ~SegmentDH();
	private:
		double dhParamA, dhParamD;
		double dhParamAlpha, dhParamTheta;
	};

	class ChainDH
	{
	public:
		/**
		 * The constructor of a chain, a new chain is always empty.
		 *
		 */
		ChainDH(size_t numSegments=0);
		ChainDH(const ChainDH& in);
		ChainDH& operator=(const ChainDH& arg);

		/**
		 * Adds a new segment to the <strong>end</strong> of the chain.
		 *
		 * @param segment The segment to add
		 */
		void addSegment(const SegmentDH& segment);

		/**
		 * Request the total number of joints in the chain.\n
		 * <strong> Important:</strong> It is not the
		 * same as the total number of segments since a segment does not
		 * need to have a joint. This function is important when
		 * creating a KDL::JntArray to use with this chain.
		 * @return total nr of joints
		 */
		unsigned int getNrOfJoints() const;
		/**
		 * Request the total number of segments in the chain.
		 * @return total number of segments
		 */
		unsigned int getNrOfSegments() const;

		/**
		 * Request the nr'd segment of the chain. There is no boundary
		 * checking.
		 *
		 * @param nr the nr of the segment starting from 0
		 *
		 * @return a constant reference to the nr'd segment
		 */
		const SegmentDH& getSegment(unsigned int nr) const;

		/**
		 * Request the nr'd segment of the chain. There is no boundary
		 * checking.
		 *
		 * @param nr the nr of the segment starting from 0
		 *
		 * @return a reference to the nr'd segment
		 */
		SegmentDH& getSegment(unsigned int nr);

		Chain getGenericChain() const;

		virtual ~ChainDH();
	private:
		unsigned int nrOfJoints;
		unsigned int nrOfSegments;
		std::vector<SegmentDH> segments;
	};

	class KinematicObservation
	{
	public:
		explicit KinematicObservation(const JntArray&, const Frame&);

		const Frame& GetToolFrame() const;
		const JntArray& GetConfiguration() const;
	private:
		Frame toolFrame;
		JntArray configuration;
	};

	class ChainParamOptimizer : public SolverI
	{
	public:
		explicit ChainParamOptimizer(const ChainDH& chain);

		void AddObservation(const KinematicObservation&);
		ChainDH OptimizeChain();
		virtual void updateInternalDataStructures();

		virtual ~ChainParamOptimizer();
	private:
		// pointer-to-implementation idiom
		struct impl;
		std::unique_ptr<impl> pImpl;
	};
}

#endif
