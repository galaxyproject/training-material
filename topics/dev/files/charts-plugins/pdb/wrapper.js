define( [ 'visualizations/myviz/pdb/plugin' ], function( pv ) {
    return Backbone.Model.extend({
        initialize: function( options ) {
            var settings = options.chart.settings;
            var viewer = pv.Viewer( document.getElementById( options.targets[ 0 ] ), {
                width       : 'auto',
                height      : 'auto',
                antialias   : true,
                outline     : true
            });
            $.ajax( {
                url     : options.dataset.download_url,
                success : function( response ) {
                    var structure = pv.io.pdb( response );
                    viewer.clear();
                    window.console.log( settings.attributes );
                    viewer.renderAs( 'protein', structure, settings.get( 'mode' ), settings.attributes );
                    viewer.centerOn( structure );
                    viewer.autoZoom();
                    options.process.resolve();
                }
            });
        }
    });
});