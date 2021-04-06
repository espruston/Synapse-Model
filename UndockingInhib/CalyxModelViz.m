classdef CalyxModelViz < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure            matlab.ui.Figure
        GridLayout          matlab.ui.container.GridLayout
        LeftPanel           matlab.ui.container.Panel
        XaxisDropDownLabel  matlab.ui.control.Label
        XaxisDropDown       matlab.ui.control.DropDown
        YaxisDropDownLabel  matlab.ui.control.Label
        YaxisDropDown       matlab.ui.control.DropDown
        ZaxisDropDownLabel  matlab.ui.control.Label
        ZaxisDropDown       matlab.ui.control.DropDown
        CenterPanel         matlab.ui.container.Panel
        UIAxes              matlab.ui.control.UIAxes
        RightPanel          matlab.ui.container.Panel
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
        twoPanelWidth = 768;
    end

    % Callbacks that handle component events
    methods (Access = private)
        
        
        function axisDropDownValueChanged(app, event)
            
            dataMat = load('CalyxBuildSet.mat').dataMat;
            
            switch app.XaxisDropDown.Value
                case 'p_{release_{docked}}'
                    xData = dataMat(:,1);
                    app.UIAxes.XLim = [0 1];
                case 'p_{release_{tethered}}'
                    xData = dataMat(:,2);
                    app.UIAxes.XLim = [0 1];
                case 'k_{docking}'
                    xData = dataMat(:,3);
                    app.UIAxes.XLim = [0 0.01];
                case 'k_{undocking}'
                    xData = dataMat(:,4);
                    app.UIAxes.XLim = [0 0.01];
                case 'k_{tether}'
                    xData = dataMat(:,5);
                    app.UIAxes.XLim = [0 0.01];
                case 'k_{untether}'
                    xData = dataMat(:,6);
                    app.UIAxes.XLim = [0 0.01];
                case 'reserve size'
                    xData = dataMat(:,7);
                    app.UIAxes.XLim = [0 100];
                case 'k_{refill}'
                    xData = dataMat(:,8);
                    app.UIAxes.XLim = [0 0.01];
                case 'C_3'
                    xData = dataMat(:,9);
                    app.UIAxes.XLim = [0 1];
                case 'C_7'
                    xData = dataMat(:,10);
                    app.UIAxes.XLim = [0 1];
            end
            
            switch app.YaxisDropDown.Value
                case 'p_{release_{docked}}'
                    yData = dataMat(:,1);
                    app.UIAxes.YLim = [0 1];
                case 'p_{release_{tethered}}'
                    yData = dataMat(:,2);
                    app.UIAxes.YLim = [0 1];
                case 'k_{docking}'
                    yData = dataMat(:,3);
                    app.UIAxes.YLim = [0 0.01];
                case 'k_{undocking}'
                    yData = dataMat(:,4);
                    app.UIAxes.YLim = [0 0.01];
                case 'k_{tether}'
                    yData = dataMat(:,5);
                    app.UIAxes.YLim = [0 0.01];
                case 'k_{untether}'
                    yData = dataMat(:,6);
                    app.UIAxes.YLim = [0 0.01];
                case 'reserve size'
                    yData = dataMat(:,7);
                    app.UIAxes.YLim = [0 100];
                case 'k_{refill}'
                    yData = dataMat(:,8);
                    app.UIAxes.YLim = [0 0.01];
                case 'C_3'
                    yData = dataMat(:,9);
                    app.UIAxes.YLim = [0 1];
                case 'C_7'
                    yData = dataMat(:,10);
                    app.UIAxes.YLim = [0 1];
            end
            
            switch app.ZaxisDropDown.Value
                case 'p_{release_{docked}}'
                    zData = dataMat(:,1);
                    app.UIAxes.ZLim = [0 1];
                case 'p_{release_{tethered}}'
                    zData = dataMat(:,2);
                    app.UIAxes.ZLim = [0 1];
                case 'k_{docking}'
                    zData = dataMat(:,3);
                    app.UIAxes.ZLim = [0 0.01];
                case 'k_{undocking}'
                    zData = dataMat(:,4);
                    app.UIAxes.ZLim = [0 0.01];
                case 'k_{tether}'
                    zData = dataMat(:,5);
                    app.UIAxes.ZLim = [0 0.01];
                case 'k_{untether}'
                    zData = dataMat(:,6);
                    app.UIAxes.ZLim = [0 0.01];
                case 'reserve size'
                    zData = dataMat(:,7);
                    app.UIAxes.ZLim = [0 100];
                case 'k_{refill}'
                    zData = dataMat(:,8);
                    app.UIAxes.ZLim = [0 0.01];
                case 'C_3'
                    zData = dataMat(:,9);
                    app.UIAxes.ZLim = [0 1];
                case 'C_7'
                    zData = dataMat(:,10);
                    app.UIAxes.ZLim = [0 1];
            end

            app.UIAxes.CLim = [min(dataMat(:,11)) max(dataMat(:,11))];
            scatter3(app.UIAxes,xData,yData,zData, 50, dataMat(:,11), 'filled', 'MarkerEdgeColor', 'k');
            colorbar(app.UIAxes,'Direction','reverse','Location','northoutside');
            
            app.UIAxes.View = [40,45];
            app.UIAxes.XGrid = 'on';
            app.UIAxes.YGrid = 'on';
            app.UIAxes.ZGrid = 'on';
            
            app.UIAxes.XDir = 'reverse';
            
            xlabel(app.UIAxes, app.XaxisDropDown.Value)
            ylabel(app.UIAxes, app.YaxisDropDown.Value)
            zlabel(app.UIAxes, app.ZaxisDropDown.Value)
            
            
            
        end
        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 3x1 grid
                app.GridLayout.RowHeight = {480, 480, 480};
                app.GridLayout.ColumnWidth = {'1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 1;
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 3;
                app.RightPanel.Layout.Column = 1;
            elseif (currentFigureWidth > app.onePanelWidth && currentFigureWidth <= app.twoPanelWidth)
                % Change to a 2x2 grid
                app.GridLayout.RowHeight = {480, 480};
                app.GridLayout.ColumnWidth = {'1x', '1x'};
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = [1,2];
                app.LeftPanel.Layout.Row = 2;
                app.LeftPanel.Layout.Column = 1;
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 2;
            else
                % Change to a 1x3 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {220, '1x', 220};
                app.LeftPanel.Layout.Row = 1;
                app.LeftPanel.Layout.Column = 1;
                app.CenterPanel.Layout.Row = 1;
                app.CenterPanel.Layout.Column = 2;
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 3;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 860 480];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {220, '1x', 220};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create XaxisDropDownLabel
            app.XaxisDropDownLabel = uilabel(app.LeftPanel);
            app.XaxisDropDownLabel.HorizontalAlignment = 'right';
            app.XaxisDropDownLabel.Position = [46 446 39 22];
            app.XaxisDropDownLabel.Text = 'X-axis';

            % Create XaxisDropDown
            app.XaxisDropDown = uidropdown(app.LeftPanel);
            app.XaxisDropDown.Items = {'p_release_docked', 'p_release_tethered', 'k_docking', 'k_undocking', 'k_tether', 'k_untether', 'reserve_size', 'k_refill', 'C_3', 'C_7'};
            app.XaxisDropDown.ItemsData = {'p_{release_{docked}}', 'p_{release_{tethered}}', 'k_{docking}', 'k_{undocking}', 'k_{tether}', 'k_{untether}', 'reserve size', 'k_{refill}', 'C_3', 'C_7'};
            app.XaxisDropDown.ValueChangedFcn = createCallbackFcn(app, @axisDropDownValueChanged, true);
            app.XaxisDropDown.Position = [100 446 100 22];
            app.XaxisDropDown.Value = 'p_{release_{docked}}';

            % Create YaxisDropDownLabel
            app.YaxisDropDownLabel = uilabel(app.LeftPanel);
            app.YaxisDropDownLabel.HorizontalAlignment = 'right';
            app.YaxisDropDownLabel.Position = [47 409 38 22];
            app.YaxisDropDownLabel.Text = 'Y-axis';

            % Create YaxisDropDown
            app.YaxisDropDown = uidropdown(app.LeftPanel);
            app.YaxisDropDown.Items = {'p_release_docked', 'p_release_tethered', 'k_docking', 'k_undocking', 'k_tether', 'k_untether', 'reserve_size', 'k_refill', 'C_3', 'C_7'};
            app.YaxisDropDown.ItemsData = {'p_{release_{docked}}', 'p_{release_{tethered}}', 'k_{docking}', 'k_{undocking}', 'k_{tether}', 'k_{untether}', 'reserve size', 'k_{refill}', 'C_3', 'C_7'};
            app.YaxisDropDown.ValueChangedFcn = createCallbackFcn(app, @axisDropDownValueChanged, true);
            app.YaxisDropDown.Position = [100 409 100 22];
            app.YaxisDropDown.Value = 'p_{release_{tethered}}';

            % Create ZaxisDropDownLabel
            app.ZaxisDropDownLabel = uilabel(app.LeftPanel);
            app.ZaxisDropDownLabel.HorizontalAlignment = 'right';
            app.ZaxisDropDownLabel.Position = [47 372 38 22];
            app.ZaxisDropDownLabel.Text = 'Z-axis';

            % Create ZaxisDropDown
            app.ZaxisDropDown = uidropdown(app.LeftPanel);
            app.ZaxisDropDown.Items = {'p_release_docked', 'p_release_tethered', 'k_docking', 'k_undocking', 'k_tether', 'k_untether', 'reserve_size', 'k_refill', 'C_3', 'C_7'};
            app.ZaxisDropDown.ItemsData = {'p_{release_{docked}}', 'p_{release_{tethered}}', 'k_{docking}', 'k_{undocking}', 'k_{tether}', 'k_{untether}', 'reserve size', 'k_{refill}', 'C_3', 'C_7'};
            app.ZaxisDropDown.ValueChangedFcn = createCallbackFcn(app, @axisDropDownValueChanged, true);
            app.ZaxisDropDown.Position = [100 372 100 22];
            app.ZaxisDropDown.Value = 'k_{docking}';

            % Create CenterPanel
            app.CenterPanel = uipanel(app.GridLayout);
            app.CenterPanel.Layout.Row = 1;
            app.CenterPanel.Layout.Column = 2;

            % Create UIAxes
            app.UIAxes = uiaxes(app.CenterPanel);
            rotate3d(app.UIAxes,'on');
            title(app.UIAxes, 'Cost vs Parameters')
            xlabel(app.UIAxes, app.XaxisDropDown.Value)
            ylabel(app.UIAxes, app.YaxisDropDown.Value)
            zlabel(app.UIAxes, app.ZaxisDropDown.Value)
            app.UIAxes.Position = [6 162 403 248];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 3;

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = CalyxModelViz

            % Create UIFigure and components
            createComponents(app)
            
            %plot defaults
            axisDropDownValueChanged(app)
            
            % Register the app with App Designer
            registerApp(app, app.UIFigure)

            if nargout == 0
                clear app
            end
        end

        % Code that executes before app deletion
        function delete(app)

            % Delete UIFigure when app is deleted
            delete(app.UIFigure)
        end
    end
end