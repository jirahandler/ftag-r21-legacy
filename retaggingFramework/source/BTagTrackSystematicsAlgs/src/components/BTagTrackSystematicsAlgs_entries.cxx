#include "GaudiKernel/DeclareFactoryEntries.h"
#include "../TrackSystematicsAlg.h"

DECLARE_NAMESPACE_ALGORITHM_FACTORY( InDet, TrackSystematicsAlg )

/** factory entries need to have the name of the package */
DECLARE_FACTORY_ENTRIES( BTagTrackSystematicsAlgs )
{
  DECLARE_NAMESPACE_ALGORITHM( InDet, TrackSystematicsAlg );
}
