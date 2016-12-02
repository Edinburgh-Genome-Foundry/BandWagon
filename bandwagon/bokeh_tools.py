from bokeh.core.properties import Dict, Int, String
from bokeh.models import TickFormatter
from bokeh.util.compiler import CoffeeScript

class FixedTickFormatter(TickFormatter):
    """
    Class used to allow custom axis tick labels on a bokeh chart
    Extends bokeh.model.formatters.TickFormatter
    """

    COFFEESCRIPT = """
        import {Model} from "model"
        import * as p from "core/properties"
        export class FixedTickFormatter extends Model
          type: 'FixedTickFormatter'
          doFormat: (ticks) ->
            labels = @get("labels")
            return (labels[tick] ? "" for tick in ticks)
          @define {
            labels: [ p.Any ]
          }
    """

    labels = Dict(Int, String, help="""
    A mapping of integer ticks values to their labels.
    """)

    __implementation__ = CoffeeScript(COFFEESCRIPT)
