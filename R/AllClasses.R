setClass("CellPopulation", representation(flow.frame="flowFrame",
                                          proportion="numeric",
                                          cell.count="numeric",
                                          channels="character", 
                                          position="logical",
                                          gates="numeric",
					  filter="matrix",
					  index="numeric"

)
)
