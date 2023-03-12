#include <exception>
#include <map>
#include <cassert>
#include "chainparamoptimizer.hpp"
#include "chainfksolverpos_recursive.hpp"
#include "chainjnttojacsolver.hpp"

namespace KDL
{
    static const Joint& dhCompatibleJointOrError(const Joint& joint)
	{
		switch (joint.getType()) {
		    case Joint::RotZ:
		    case Joint::TransZ:
		    case Joint::Fixed:
			    return joint;
		    default:
			    throw std::invalid_argument("DH convention requires joint of type Joint::RotZ or Joint::TransZ");
		}
	}

	SegmentDH::SegmentDH(const std::string& name, const Joint& joint, double a, double alpha, double d, double theta, const RigidBodyInertia& I)
	      : Segment(name, dhCompatibleJointOrError(joint), Frame::DH_Craig1989(a, alpha, d, theta), I),
	        dhParamA(a), dhParamD(d), dhParamAlpha(), dhParamTheta(theta)
	{

	}

	SegmentDH::SegmentDH(const SegmentDH &other)
	    : Segment(other),
	      dhParamA(other.dhParamA), dhParamD(other.dhParamA),
	      dhParamAlpha(other.dhParamAlpha), dhParamTheta(other.dhParamTheta)
	{
	}

	SegmentDH& SegmentDH::operator =(const SegmentDH& arg)
	{
		dhParamA = arg.dhParamA;
		dhParamD = arg.dhParamD;
		dhParamAlpha = arg.dhParamAlpha;
		dhParamTheta = arg.dhParamTheta;
		Segment::operator =(arg);
		return *this;
	}

	SegmentDH::~SegmentDH()
	{
	}

	ChainDH::ChainDH(size_t numSegments)
	{
		if (segments.capacity()<numSegments)
			segments.reserve(numSegments);
	}

	ChainDH::ChainDH(const ChainDH &in)
	    : segments(in.segments)
	{
		nrOfJoints = in.nrOfJoints;
		nrOfSegments = in.nrOfSegments;
	}

	ChainDH::~ChainDH()
	{
	}

	ChainDH& ChainDH::operator = (const ChainDH& rhs)
	{
		nrOfJoints = rhs.nrOfJoints;
		nrOfSegments = rhs.nrOfSegments;
		segments.clear();
		for (const SegmentDH& seg: rhs.segments)
		{
			this->segments.push_back(seg);
		}
		return *this;
	}

	void ChainDH::addSegment(const SegmentDH &segment)
	{
		if(Joint::Fixed==segment.getJoint().getType() && (0!=nrOfSegments))
		{
			SegmentDH lastSeg = segments.back();

			if (lastSeg.getJoint().getType() == Joint::Fixed)
				throw std::invalid_argument("DH chain expects one joint per segement except for base and tool links");
		}

		segments.push_back(segment);
		nrOfSegments++;

		if(Joint::Fixed!=segment.getJoint().getType())
			nrOfJoints++;
	}

	unsigned int ChainDH::getNrOfJoints() const
	{
		return nrOfJoints;
	}

	unsigned int ChainDH::getNrOfSegments() const
	{
		return nrOfSegments;
	}

	const SegmentDH& ChainDH::getSegment(unsigned int nr) const
	{
		return segments[nr];
	}

	Chain ChainDH::getGenericChain() const
	{
		Chain c;
		for (const SegmentDH& seg : segments) {
			c.addSegment(seg);
		}
		return c;
	}

	SegmentDH& ChainDH::getSegment(unsigned int nr)
	{
		return segments[nr];
	}


	KinematicObservation::KinematicObservation(const JntArray & config, const Frame & tool)
	    : toolFrame(tool),
	      configuration(config)
	{
	}

	const JntArray& KinematicObservation::GetConfiguration() const
	{
		return configuration;
	}

	const Frame& KinematicObservation::GetToolFrame() const
	{
		return toolFrame;
	}

	struct ChainParamOptimizer::impl
	{
		ChainDH initialChain;
		std::vector<KinematicObservation> observations;
		ChainDH optimizedChain;

		impl(const ChainDH & chain);
		ChainDH OptimizeChain();
		void updateInternalDataStructures();
	};

	struct SegmentOptimizationMetadata
	{
		int idxOfFirstParam;
		int numParams;
	};

	std::vector<SegmentOptimizationMetadata> build_metadata(const ChainDH& chain)
	{
		const int nrOfSegments = chain.getNrOfSegments();
		std::vector<SegmentOptimizationMetadata> chainMeta(nrOfSegments);

		size_t paramCounter = 0;
		for (int segIdx = 0; segIdx < nrOfSegments; ++segIdx) {
			auto seg = chain.getSegment(segIdx);

			SegmentOptimizationMetadata segMeta;
			segMeta.numParams = 4; // at the moment: optimize 4 DH params
			segMeta.idxOfFirstParam = paramCounter;
			paramCounter += segMeta.numParams;

			chainMeta[segIdx] = segMeta;
		}

		return chainMeta;
	}

	ChainParamOptimizer::ChainParamOptimizer(const ChainDH & chain)
	    : pImpl(std::make_unique<ChainParamOptimizer::impl>(chain))
	{
	}

	ChainParamOptimizer::~ChainParamOptimizer()
	{
	}

	void ChainParamOptimizer::AddObservation(const KinematicObservation & newObservation)
	{
		pImpl->observations.push_back(newObservation);
	}

	ChainDH ChainParamOptimizer::OptimizeChain()
	{
		return pImpl->OptimizeChain();
	}

	void ChainParamOptimizer::updateInternalDataStructures() {
		pImpl->updateInternalDataStructures();
	}

	ChainParamOptimizer::impl::impl(const ChainDH &chain)
	    : initialChain(chain)
	{
	}

	ChainDH ChainParamOptimizer::impl::OptimizeChain()
	{
		if (0>=initialChain.getNrOfSegments())
			throw std::runtime_error("chain does not contain any segmetns");
		if (0>=initialChain.getNrOfJoints())
			throw std::runtime_error("chain does not contain any joints");

		auto segMeta = build_metadata(initialChain);
		auto fkPosSolver = ChainFkSolverPos_recursive(initialChain.getGenericChain());

		assert(0 < segMeta.size());

		const int numParams = segMeta.back().idxOfFirstParam + segMeta.back().numParams;
		const int numObservations = observations.size();
		const int numConstraints = numObservations * 6;

		if (0 >= numParams) {
			throw new std::runtime_error("number of optimization parameters is zero");
			// Todo: throw exception or just return a copy of the input chain?
		}
		if (numConstraints < numParams) {
			throw new std::runtime_error("too few observations: underdetermined system of equations!");
		}

		optimizedChain = ChainDH(initialChain);
		return optimizedChain;
	}

	void ChainParamOptimizer::impl::updateInternalDataStructures()
	{
ex
	}
}
