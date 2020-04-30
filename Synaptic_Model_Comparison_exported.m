classdef Synaptic_Model_Comparison_exported < matlab.apps.AppBase

    % Properties that correspond to app components
    properties (Access = public)
        UIFigure                        matlab.ui.Figure
        GridLayout                      matlab.ui.container.GridLayout
        LeftPanel                       matlab.ui.container.Panel
        TabGroup                        matlab.ui.container.TabGroup
        VesicleMaturationTab            matlab.ui.container.Tab
        FrequencyHzEditFieldLabel       matlab.ui.control.Label
        FrequencyHzEditField            matlab.ui.control.NumericEditField
        Stage01CharacteristicTimesEditFieldLabel  matlab.ui.control.Label
        Stage01CharacteristicTimesEditField  matlab.ui.control.NumericEditField
        Stage12CharacteristicTimesEditFieldLabel  matlab.ui.control.Label
        Stage12CharacteristicTimesEditField  matlab.ui.control.NumericEditField
        Stage23CharacteristicTimesEditFieldLabel  matlab.ui.control.Label
        Stage23CharacteristicTimesEditField  matlab.ui.control.NumericEditField
        ProbabilityofReleaseStage1EditFieldLabel  matlab.ui.control.Label
        ProbabilityofReleaseStage1EditField  matlab.ui.control.NumericEditField
        ProbabilityofReleaseStage2EditFieldLabel  matlab.ui.control.Label
        ProbabilityofReleaseStage2EditField  matlab.ui.control.NumericEditField
        ProbabilityofReleaseStage3EditFieldLabel  matlab.ui.control.Label
        ProbabilityofReleaseStage3EditField  matlab.ui.control.NumericEditField
        NumberofVesiclestotalEditFieldLabel  matlab.ui.control.Label
        NumberofVesiclestotalEditField  matlab.ui.control.NumericEditField
        NumberofPulsesEditFieldLabel    matlab.ui.control.Label
        NumberofPulsesEditField         matlab.ui.control.NumericEditField
        CurrentperVesicleEditFieldLabel  matlab.ui.control.Label
        CurrentperVesicleEditField      matlab.ui.control.NumericEditField
        LooselyFillinganemptydockingsiteLabel  matlab.ui.control.Label
        Stage1becomingstage2Label       matlab.ui.control.Label
        Stage2becomingfacilitatedLabel  matlab.ui.control.Label
        DependentVariableDropDownLabel  matlab.ui.control.Label
        DependentVariableDropDown       matlab.ui.control.DropDown
        SetParametersButton             matlab.ui.control.Button
        SetparametersbeforegraphingLabel  matlab.ui.control.Label
        RefreshButton                   matlab.ui.control.Button
        UITable                         matlab.ui.control.Table
        TwoPoolTab                      matlab.ui.container.Tab
        FrequencyHzEditField_2Label     matlab.ui.control.Label
        FrequencyHzEditField_2          matlab.ui.control.NumericEditField
        ProbabilityofReleasePool1EditFieldLabel  matlab.ui.control.Label
        ProbabilityofReleasePool1EditField  matlab.ui.control.NumericEditField
        ProbabilityofReleasePool2Label  matlab.ui.control.Label
        ProbabilityofReleasePool2EditField  matlab.ui.control.NumericEditField
        ProbabilityofReleaseFacilitatedEditFieldLabel  matlab.ui.control.Label
        ProbabilityofReleaseFacilitatedEditField  matlab.ui.control.NumericEditField
        NumberofVesiclesPool1EditFieldLabel  matlab.ui.control.Label
        NumberofVesiclesPool1EditField  matlab.ui.control.NumericEditField
        NumberofPulsesEditField_2Label  matlab.ui.control.Label
        NumberofPulsesEditField_2       matlab.ui.control.NumericEditField
        CurrentperVesicleEditField_2Label  matlab.ui.control.Label
        CurrentperVesicleEditField_2    matlab.ui.control.NumericEditField
        FacilitationCharacteristicTimesEditFieldLabel  matlab.ui.control.Label
        FacilitationCharacteristicTimesEditField  matlab.ui.control.NumericEditField
        Pool2RefillCharacteristicTimesEditFieldLabel  matlab.ui.control.Label
        Pool2RefillCharacteristicTimesEditField  matlab.ui.control.NumericEditField
        Pool1RefillCharacteristicTimesEditFieldLabel  matlab.ui.control.Label
        Pool1RefillCharacteristicTimesEditField  matlab.ui.control.NumericEditField
        NumberofVesiclesPool2EditFieldLabel  matlab.ui.control.Label
        NumberofVesiclesPool2EditField  matlab.ui.control.NumericEditField
        RapidResponsePoolLabel          matlab.ui.control.Label
        SteadyStatePoolLabel            matlab.ui.control.Label
        SetParametersButton_2           matlab.ui.control.Button
        SetparametersbeforegraphingLabel_2  matlab.ui.control.Label
        RefreshButton_2                 matlab.ui.control.Button
        DependentVariableDropDown_2Label  matlab.ui.control.Label
        DependentVariableDropDown_2     matlab.ui.control.DropDown
        UITable2                        matlab.ui.control.Table
        RightPanel                      matlab.ui.container.Panel
        UIAxes                          matlab.ui.control.UIAxes
    end

    % Properties that correspond to apps with auto-reflow
    properties (Access = private)
        onePanelWidth = 576;
    end

    % Callbacks that handle component events
    methods (Access = private)

        % Value changed function: FrequencyHzEditField
        function FrequencyHzEditFieldValueChanged(app, event)

            
        end

        % Button pushed function: SetParametersButton
        function SetParametersButtonPushed(app, event)
            app.UITable.Data = Vesicle_maturation_model_func(app.FrequencyHzEditField.Value, app.Stage01CharacteristicTimesEditField.Value, app.Stage12CharacteristicTimesEditField.Value, app.Stage23CharacteristicTimesEditField.Value, app.ProbabilityofReleaseStage1EditField.Value, app.ProbabilityofReleaseStage2EditField.Value, app.ProbabilityofReleaseStage3EditField.Value, app.NumberofVesiclestotalEditField.Value, app.NumberofPulsesEditField.Value, app.CurrentperVesicleEditField.Value);    
        end

        % Value changed function: DependentVariableDropDown
        function DependentVariableDropDownValueChanged(app, event)
            value = app.DependentVariableDropDown.Value;
            value = str2double(value);
            %f = msgbox(sprintf('%s, %c', class(value), value));
            title_var = ["Empty Docking Sites", "Stage 1 Vesicles", "Stage 2 Vesicles", "Stage 3 Vesicles", "Stage 1 Vesicles Released", "Stage 2 Vesicles Released", "Stage 3 Vesicles Released", "Total Vesicles Released", "Current", "New Docking Sites", "New Stage 1 Vesicles", "New Stage 2 Vesicles", "New Stage 3 Vesicles"];
            app.UIAxes.Title.String = sprintf('%s at %d Hz', title_var(value-1), app.FrequencyHzEditField.Value);
            app.UIAxes.XLabel.String = 'Pulse Number';
            app.UIAxes.XLim = [0, app.NumberofPulsesEditField.Value];
            app.UIAxes.YLabel.String = sprintf('%s', title_var(value-1));
            %app.UIAxes.YLim = [min(app.UITable.Data(:,value)), max(app.UITable.Data(:,value))];
            bar(app.UIAxes, app.UITable.Data(:,1), app.UITable.Data(:,value));
        end

        % Button pushed function: RefreshButton
        function RefreshButtonPushed(app, event)
            DependentVariableDropDownValueChanged(app, event);
        end

        % Button pushed function: SetParametersButton_2
        function SetParametersButton_2Pushed(app, event)
            app.UITable2.Data = Two_Pool_Model_func(app.FrequencyHzEditField_2.Value, app.ProbabilityofReleasePool1EditField.Value, app.ProbabilityofReleasePool2EditField.Value, app.ProbabilityofReleaseFacilitatedEditField.Value, app.ProbabilityofReleaseFacilitatedEditField.Value, app.NumberofVesiclesPool1EditField.Value, app.NumberofVesiclesPool2EditField.Value, app.NumberofPulsesEditField_2.Value, app.CurrentperVesicleEditField_2.Value, app.Pool1RefillCharacteristicTimesEditField.Value, app.Pool2RefillCharacteristicTimesEditField.Value, app.FacilitationCharacteristicTimesEditField.Value);
        end

        % Value changed function: DependentVariableDropDown_2
        function DependentVariableDropDown_2ValueChanged(app, event)
            value = app.DependentVariableDropDown_2.Value;
            value = str2double(value);
            %f = msgbox(sprintf('%s, %c', class(value), value));
            title_var = ["Empty Docking Sites (Pool 1)", "Empty Docking Sites (Pool 2)", "Pool 1 Vesicles", "Pool 2 Vesicles", "Facilitated Pool 1 Vesicles", "Facilitated Pool 2 Vesicles", "Vesicles Released From Pool 1", "Vesicles Released From Pool 2", "Facilitated Vesicles Released From Pool 1", "Facilitated Vesicles Released From Pool 2", "Total Vesicles Released", "Current", "New Pool 1 Vesicles", "New Pool 2 Vesicles", "New Facilitated Pool 1 Vesicles", "New Facilitated Pool 2 Vesicles"];
            app.UIAxes.Title.String = sprintf('%s at %d Hz', title_var(value-1), app.FrequencyHzEditField_2.Value);
            app.UIAxes.XLabel.String = 'Pulse Number';
            app.UIAxes.XLim = [0, app.NumberofPulsesEditField_2.Value];
            app.UIAxes.YLabel.String = sprintf('%s', title_var(value-1));
            %app.UIAxes.YLim = [min(app.UITable2.Data(:,value)), max(app.UITable2.Data(:,value))];
            bar(app.UIAxes, app.UITable2.Data(:,1), app.UITable2.Data(:,value));
        end

        % Button pushed function: RefreshButton_2
        function RefreshButton_2Pushed(app, event)
            DependentVariableDropDown_2ValueChanged(app, event);
        end

        % Changes arrangement of the app based on UIFigure width
        function updateAppLayout(app, event)
            currentFigureWidth = app.UIFigure.Position(3);
            if(currentFigureWidth <= app.onePanelWidth)
                % Change to a 2x1 grid
                app.GridLayout.RowHeight = {701, 701};
                app.GridLayout.ColumnWidth = {'1x'};
                app.RightPanel.Layout.Row = 2;
                app.RightPanel.Layout.Column = 1;
            else
                % Change to a 1x2 grid
                app.GridLayout.RowHeight = {'1x'};
                app.GridLayout.ColumnWidth = {549, '1x'};
                app.RightPanel.Layout.Row = 1;
                app.RightPanel.Layout.Column = 2;
            end
        end
    end

    % Component initialization
    methods (Access = private)

        % Create UIFigure and components
        function createComponents(app)
            
            frequency_vm=100;
            ct01_vm=0.2;
            ct12_vm=7.5;
            ct23_vm=0.025;
            p1_vm=0.028;
            p2_vm=0.09;
            p3_vm=0.09;
            nv_vm=100;
            np_vm=100;
            cpv_vm=-1;

            % Create UIFigure and hide until all components are created
            app.UIFigure = uifigure('Visible', 'off');
            app.UIFigure.AutoResizeChildren = 'off';
            app.UIFigure.Position = [100 100 1230 701];
            app.UIFigure.Name = 'MATLAB App';
            app.UIFigure.SizeChangedFcn = createCallbackFcn(app, @updateAppLayout, true);

            % Create GridLayout
            app.GridLayout = uigridlayout(app.UIFigure);
            app.GridLayout.ColumnWidth = {549, '1x'};
            app.GridLayout.RowHeight = {'1x'};
            app.GridLayout.ColumnSpacing = 0;
            app.GridLayout.RowSpacing = 0;
            app.GridLayout.Padding = [0 0 0 0];
            app.GridLayout.Scrollable = 'on';

            % Create LeftPanel
            app.LeftPanel = uipanel(app.GridLayout);
            app.LeftPanel.Layout.Row = 1;
            app.LeftPanel.Layout.Column = 1;

            % Create TabGroup
            app.TabGroup = uitabgroup(app.LeftPanel);
            app.TabGroup.Position = [6 6 538 694];

            % Create VesicleMaturationTab
            app.VesicleMaturationTab = uitab(app.TabGroup);
            app.VesicleMaturationTab.Title = 'Vesicle Maturation';

            % Create FrequencyHzEditFieldLabel
            app.FrequencyHzEditFieldLabel = uilabel(app.VesicleMaturationTab);
            app.FrequencyHzEditFieldLabel.HorizontalAlignment = 'right';
            app.FrequencyHzEditFieldLabel.Position = [1 647 88 22];
            app.FrequencyHzEditFieldLabel.Text = 'Frequency (Hz)';

            % Create FrequencyHzEditField
            app.FrequencyHzEditField = uieditfield(app.VesicleMaturationTab, 'numeric');
            app.FrequencyHzEditField.Limits = [0 Inf];
            app.FrequencyHzEditField.ValueChangedFcn = createCallbackFcn(app, @FrequencyHzEditFieldValueChanged, true);
            app.FrequencyHzEditField.Position = [212 647 100 22];
            app.FrequencyHzEditField.Value = frequency_vm;

            % Create Stage01CharacteristicTimesEditFieldLabel
            app.Stage01CharacteristicTimesEditFieldLabel = uilabel(app.VesicleMaturationTab);
            app.Stage01CharacteristicTimesEditFieldLabel.HorizontalAlignment = 'right';
            app.Stage01CharacteristicTimesEditFieldLabel.Position = [1 626 196 22];
            app.Stage01CharacteristicTimesEditFieldLabel.Text = 'Stage 0 -> 1 Characteristic Time (s)';

            % Create Stage01CharacteristicTimesEditField
            app.Stage01CharacteristicTimesEditField = uieditfield(app.VesicleMaturationTab, 'numeric');
            app.Stage01CharacteristicTimesEditField.Limits = [0 Inf];
            app.Stage01CharacteristicTimesEditField.Position = [212 626 100 22];
            app.Stage01CharacteristicTimesEditField.Value = ct01_vm;

            % Create Stage12CharacteristicTimesEditFieldLabel
            app.Stage12CharacteristicTimesEditFieldLabel = uilabel(app.VesicleMaturationTab);
            app.Stage12CharacteristicTimesEditFieldLabel.HorizontalAlignment = 'right';
            app.Stage12CharacteristicTimesEditFieldLabel.Position = [1 605 196 22];
            app.Stage12CharacteristicTimesEditFieldLabel.Text = 'Stage 1 -> 2 Characteristic Time (s)';

            % Create Stage12CharacteristicTimesEditField
            app.Stage12CharacteristicTimesEditField = uieditfield(app.VesicleMaturationTab, 'numeric');
            app.Stage12CharacteristicTimesEditField.Limits = [0 Inf];
            app.Stage12CharacteristicTimesEditField.Position = [212 605 100 22];
            app.Stage12CharacteristicTimesEditField.Value = ct12_vm;

            % Create Stage23CharacteristicTimesEditFieldLabel
            app.Stage23CharacteristicTimesEditFieldLabel = uilabel(app.VesicleMaturationTab);
            app.Stage23CharacteristicTimesEditFieldLabel.HorizontalAlignment = 'right';
            app.Stage23CharacteristicTimesEditFieldLabel.Position = [0 584 196 22];
            app.Stage23CharacteristicTimesEditFieldLabel.Text = 'Stage 2 -> 3 Characteristic Time (s)';

            % Create Stage23CharacteristicTimesEditField
            app.Stage23CharacteristicTimesEditField = uieditfield(app.VesicleMaturationTab, 'numeric');
            app.Stage23CharacteristicTimesEditField.Limits = [0 Inf];
            app.Stage23CharacteristicTimesEditField.Position = [212 584 100 22];
            app.Stage23CharacteristicTimesEditField.Value = ct23_vm;

            % Create ProbabilityofReleaseStage1EditFieldLabel
            app.ProbabilityofReleaseStage1EditFieldLabel = uilabel(app.VesicleMaturationTab);
            app.ProbabilityofReleaseStage1EditFieldLabel.HorizontalAlignment = 'right';
            app.ProbabilityofReleaseStage1EditFieldLabel.Position = [0 563 175 22];
            app.ProbabilityofReleaseStage1EditFieldLabel.Text = 'Probability of Release (Stage 1)';

            % Create ProbabilityofReleaseStage1EditField
            app.ProbabilityofReleaseStage1EditField = uieditfield(app.VesicleMaturationTab, 'numeric');
            app.ProbabilityofReleaseStage1EditField.Limits = [0 1];
            app.ProbabilityofReleaseStage1EditField.Position = [212 563 100 22];
            app.ProbabilityofReleaseStage1EditField.Value = p1_vm;

            % Create ProbabilityofReleaseStage2EditFieldLabel
            app.ProbabilityofReleaseStage2EditFieldLabel = uilabel(app.VesicleMaturationTab);
            app.ProbabilityofReleaseStage2EditFieldLabel.HorizontalAlignment = 'right';
            app.ProbabilityofReleaseStage2EditFieldLabel.Position = [0 542 175 22];
            app.ProbabilityofReleaseStage2EditFieldLabel.Text = 'Probability of Release (Stage 2)';

            % Create ProbabilityofReleaseStage2EditField
            app.ProbabilityofReleaseStage2EditField = uieditfield(app.VesicleMaturationTab, 'numeric');
            app.ProbabilityofReleaseStage2EditField.Limits = [0 1];
            app.ProbabilityofReleaseStage2EditField.Position = [212 542 100 22];
            app.ProbabilityofReleaseStage2EditField.Value = p2_vm;

            % Create ProbabilityofReleaseStage3EditFieldLabel
            app.ProbabilityofReleaseStage3EditFieldLabel = uilabel(app.VesicleMaturationTab);
            app.ProbabilityofReleaseStage3EditFieldLabel.HorizontalAlignment = 'right';
            app.ProbabilityofReleaseStage3EditFieldLabel.Position = [1 521 175 22];
            app.ProbabilityofReleaseStage3EditFieldLabel.Text = 'Probability of Release (Stage 3)';

            % Create ProbabilityofReleaseStage3EditField
            app.ProbabilityofReleaseStage3EditField = uieditfield(app.VesicleMaturationTab, 'numeric');
            app.ProbabilityofReleaseStage3EditField.Limits = [0 1];
            app.ProbabilityofReleaseStage3EditField.Position = [212 521 100 22];
            app.ProbabilityofReleaseStage3EditField.Value = p3_vm;

            % Create NumberofVesiclestotalEditFieldLabel
            app.NumberofVesiclestotalEditFieldLabel = uilabel(app.VesicleMaturationTab);
            app.NumberofVesiclestotalEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofVesiclestotalEditFieldLabel.Position = [1 500 143 22];
            app.NumberofVesiclestotalEditFieldLabel.Text = 'Number of Vesicles (total)';

            % Create NumberofVesiclestotalEditField
            app.NumberofVesiclestotalEditField = uieditfield(app.VesicleMaturationTab, 'numeric');
            app.NumberofVesiclestotalEditField.Limits = [0 Inf];
            app.NumberofVesiclestotalEditField.Position = [212 500 100 22];
            app.NumberofVesiclestotalEditField.Value = nv_vm;

            % Create NumberofPulsesEditFieldLabel
            app.NumberofPulsesEditFieldLabel = uilabel(app.VesicleMaturationTab);
            app.NumberofPulsesEditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofPulsesEditFieldLabel.Position = [1 479 101 22];
            app.NumberofPulsesEditFieldLabel.Text = 'Number of Pulses';

            % Create NumberofPulsesEditField
            app.NumberofPulsesEditField = uieditfield(app.VesicleMaturationTab, 'numeric');
            app.NumberofPulsesEditField.Limits = [0 Inf];
            app.NumberofPulsesEditField.Position = [212 479 100 22];
            app.NumberofPulsesEditField.Value = np_vm;

            % Create CurrentperVesicleEditFieldLabel
            app.CurrentperVesicleEditFieldLabel = uilabel(app.VesicleMaturationTab);
            app.CurrentperVesicleEditFieldLabel.HorizontalAlignment = 'right';
            app.CurrentperVesicleEditFieldLabel.Position = [1 458 108 22];
            app.CurrentperVesicleEditFieldLabel.Text = 'Current per Vesicle';

            % Create CurrentperVesicleEditField
            app.CurrentperVesicleEditField = uieditfield(app.VesicleMaturationTab, 'numeric');
            app.CurrentperVesicleEditField.Position = [212 458 100 22];
            app.CurrentperVesicleEditField.Value = cpv_vm;

            % Create LooselyFillinganemptydockingsiteLabel
            app.LooselyFillinganemptydockingsiteLabel = uilabel(app.VesicleMaturationTab);
            app.LooselyFillinganemptydockingsiteLabel.Position = [317 626 210 22];
            app.LooselyFillinganemptydockingsiteLabel.Text = '"Loosely Filling" an empty docking site';

            % Create Stage1becomingstage2Label
            app.Stage1becomingstage2Label = uilabel(app.VesicleMaturationTab);
            app.Stage1becomingstage2Label.Position = [317 605 145 22];
            app.Stage1becomingstage2Label.Text = 'Stage 1 becoming stage 2';

            % Create Stage2becomingfacilitatedLabel
            app.Stage2becomingfacilitatedLabel = uilabel(app.VesicleMaturationTab);
            app.Stage2becomingfacilitatedLabel.Position = [317 584 156 22];
            app.Stage2becomingfacilitatedLabel.Text = 'Stage 2 becoming facilitated';

            % Create DependentVariableDropDownLabel
            app.DependentVariableDropDownLabel = uilabel(app.VesicleMaturationTab);
            app.DependentVariableDropDownLabel.HorizontalAlignment = 'right';
            app.DependentVariableDropDownLabel.Position = [0 374 111 22];
            app.DependentVariableDropDownLabel.Text = 'Dependent Variable';

            % Create DependentVariableDropDown
            app.DependentVariableDropDown = uidropdown(app.VesicleMaturationTab);
            app.DependentVariableDropDown.Items = {'Empty Docking Sites', 'Stage 1 Vesicles', 'Stage 2 Vesicles', 'Stage 3 Vesicles', 'Stage 1 Vesicles Released', 'Stage 2 Vesicles Released', 'Stage 3 Vesicles Released', 'Total Vesicles Released', 'Current', 'New Empty Docking Sites', 'New Stage 1 Vesicles', 'New Stage 2 Vesicles', 'New Stage 3 Vesicles'};
            app.DependentVariableDropDown.ItemsData = {'2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14'};
            app.DependentVariableDropDown.ValueChangedFcn = createCallbackFcn(app, @DependentVariableDropDownValueChanged, true);
            app.DependentVariableDropDown.Position = [126 374 224 22];
            app.DependentVariableDropDown.Value = '2';

            % Create SetParametersButton
            app.SetParametersButton = uibutton(app.VesicleMaturationTab, 'push');
            app.SetParametersButton.ButtonPushedFcn = createCallbackFcn(app, @SetParametersButtonPushed, true);
            app.SetParametersButton.Position = [6 416 100 22];
            app.SetParametersButton.Text = 'Set Parameters';

            % Create SetparametersbeforegraphingLabel
            app.SetparametersbeforegraphingLabel = uilabel(app.VesicleMaturationTab);
            app.SetparametersbeforegraphingLabel.FontWeight = 'bold';
            app.SetparametersbeforegraphingLabel.Position = [126 416 187 22];
            app.SetparametersbeforegraphingLabel.Text = 'Set parameters before graphing';

            % Create RefreshButton
            app.RefreshButton = uibutton(app.VesicleMaturationTab, 'push');
            app.RefreshButton.ButtonPushedFcn = createCallbackFcn(app, @RefreshButtonPushed, true);
            app.RefreshButton.Position = [362 374 100 22];
            app.RefreshButton.Text = 'Refresh';

            % Create UITable
            app.UITable = uitable(app.VesicleMaturationTab);
            app.UITable.ColumnName = {'Pulse Number'; 'Empty Docking Sites'; 'Stage 1 Vesicles'; 'Stage 2 Vesicles'; 'Stage 3 Vesicles'; 'Stage 1 Vesicles Released'; 'Stage 2 Vesicles Released'; 'Stage 3 Vesicles Released'; 'Total Vesicles Released'; 'Current'; 'New Docking Sites'; 'New Stage 1 Vesicles'; 'New Stage 2 Vesicles'; 'New Stage 3 Vesicles'};
            app.UITable.RowName = {};
            app.UITable.Position = [1 1 536 363];

            % Create TwoPoolTab
            app.TwoPoolTab = uitab(app.TabGroup);
            app.TwoPoolTab.Title = 'Two Pool';

            frequency_2p=100;
            ct1_2p=7.5;
            ct2_2p=0.28;
            ctf_2p=0.125;
            p1_2p=0.14;
            p2_2p=0.023;
            pf_2p=0.035;
            nv1_2p=30;
            nv2_2p=70;
            np_2p=100;
            cpv_2p=-1;

            % Create FrequencyHzEditField_2Label
            app.FrequencyHzEditField_2Label = uilabel(app.TwoPoolTab);
            app.FrequencyHzEditField_2Label.HorizontalAlignment = 'right';
            app.FrequencyHzEditField_2Label.Position = [0 647 88 22];
            app.FrequencyHzEditField_2Label.Text = 'Frequency (Hz)';

            % Create FrequencyHzEditField_2
            app.FrequencyHzEditField_2 = uieditfield(app.TwoPoolTab, 'numeric');
            app.FrequencyHzEditField_2.Limits = [0 Inf];
            app.FrequencyHzEditField_2.Position = [242 647 99 22];
            app.FrequencyHzEditField_2.Value = frequency_2p;

            % Create ProbabilityofReleasePool1EditFieldLabel
            app.ProbabilityofReleasePool1EditFieldLabel = uilabel(app.TwoPoolTab);
            app.ProbabilityofReleasePool1EditFieldLabel.HorizontalAlignment = 'right';
            app.ProbabilityofReleasePool1EditFieldLabel.Position = [0 626 168 22];
            app.ProbabilityofReleasePool1EditFieldLabel.Text = 'Probability of Release (Pool 1)';

            % Create ProbabilityofReleasePool1EditField
            app.ProbabilityofReleasePool1EditField = uieditfield(app.TwoPoolTab, 'numeric');
            app.ProbabilityofReleasePool1EditField.Limits = [0 1];
            app.ProbabilityofReleasePool1EditField.Position = [242 626 99 22];
            app.ProbabilityofReleasePool1EditField.Value = p1_2p;

            % Create ProbabilityofReleasePool2Label
            app.ProbabilityofReleasePool2Label = uilabel(app.TwoPoolTab);
            app.ProbabilityofReleasePool2Label.HorizontalAlignment = 'right';
            app.ProbabilityofReleasePool2Label.Position = [0 605 168 22];
            app.ProbabilityofReleasePool2Label.Text = 'Probability of Release (Pool 2)';

            % Create ProbabilityofReleasePool2EditField
            app.ProbabilityofReleasePool2EditField = uieditfield(app.TwoPoolTab, 'numeric');
            app.ProbabilityofReleasePool2EditField.Limits = [0 1];
            app.ProbabilityofReleasePool2EditField.Position = [242 605 99 22];
            app.ProbabilityofReleasePool2EditField.Value = p2_2p;

            % Create ProbabilityofReleaseFacilitatedEditFieldLabel
            app.ProbabilityofReleaseFacilitatedEditFieldLabel = uilabel(app.TwoPoolTab);
            app.ProbabilityofReleaseFacilitatedEditFieldLabel.HorizontalAlignment = 'right';
            app.ProbabilityofReleaseFacilitatedEditFieldLabel.Position = [0 584 188 22];
            app.ProbabilityofReleaseFacilitatedEditFieldLabel.Text = 'Probability of Release (Facilitated)';

            % Create ProbabilityofReleaseFacilitatedEditField
            app.ProbabilityofReleaseFacilitatedEditField = uieditfield(app.TwoPoolTab, 'numeric');
            app.ProbabilityofReleaseFacilitatedEditField.Limits = [0 1];
            app.ProbabilityofReleaseFacilitatedEditField.Position = [242 584 99 22];
            app.ProbabilityofReleaseFacilitatedEditField.Value = pf_2p;

            % Create NumberofVesiclesPool1EditFieldLabel
            app.NumberofVesiclesPool1EditFieldLabel = uilabel(app.TwoPoolTab);
            app.NumberofVesiclesPool1EditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofVesiclesPool1EditFieldLabel.Position = [1 563 154 22];
            app.NumberofVesiclesPool1EditFieldLabel.Text = 'Number of Vesicles (Pool 1)';

            % Create NumberofVesiclesPool1EditField
            app.NumberofVesiclesPool1EditField = uieditfield(app.TwoPoolTab, 'numeric');
            app.NumberofVesiclesPool1EditField.Limits = [0 Inf];
            app.NumberofVesiclesPool1EditField.Position = [242 563 99 22];
            app.NumberofVesiclesPool1EditField.Value = nv1_2p;

            % Create NumberofPulsesEditField_2Label
            app.NumberofPulsesEditField_2Label = uilabel(app.TwoPoolTab);
            app.NumberofPulsesEditField_2Label.HorizontalAlignment = 'right';
            app.NumberofPulsesEditField_2Label.Position = [1 521 101 22];
            app.NumberofPulsesEditField_2Label.Text = 'Number of Pulses';

            % Create NumberofPulsesEditField_2
            app.NumberofPulsesEditField_2 = uieditfield(app.TwoPoolTab, 'numeric');
            app.NumberofPulsesEditField_2.Limits = [0 Inf];
            app.NumberofPulsesEditField_2.Position = [242 521 99 22];
            app.NumberofPulsesEditField_2.Value = np_2p;

            % Create CurrentperVesicleEditField_2Label
            app.CurrentperVesicleEditField_2Label = uilabel(app.TwoPoolTab);
            app.CurrentperVesicleEditField_2Label.HorizontalAlignment = 'right';
            app.CurrentperVesicleEditField_2Label.Position = [0 500 108 22];
            app.CurrentperVesicleEditField_2Label.Text = 'Current per Vesicle';

            % Create CurrentperVesicleEditField_2
            app.CurrentperVesicleEditField_2 = uieditfield(app.TwoPoolTab, 'numeric');
            app.CurrentperVesicleEditField_2.Position = [242 500 99 22];
            app.CurrentperVesicleEditField_2.Value = cpv_2p;

            % Create FacilitationCharacteristicTimesEditFieldLabel
            app.FacilitationCharacteristicTimesEditFieldLabel = uilabel(app.TwoPoolTab);
            app.FacilitationCharacteristicTimesEditFieldLabel.HorizontalAlignment = 'right';
            app.FacilitationCharacteristicTimesEditFieldLabel.Position = [1 437 186 22];
            app.FacilitationCharacteristicTimesEditFieldLabel.Text = 'Facilitation Characteristic Time (s)';

            % Create FacilitationCharacteristicTimesEditField
            app.FacilitationCharacteristicTimesEditField = uieditfield(app.TwoPoolTab, 'numeric');
            app.FacilitationCharacteristicTimesEditField.Limits = [0 Inf];
            app.FacilitationCharacteristicTimesEditField.Position = [242 437 99 22];
            app.FacilitationCharacteristicTimesEditField.Value = ctf_2p;

            % Create Pool2RefillCharacteristicTimesEditFieldLabel
            app.Pool2RefillCharacteristicTimesEditFieldLabel = uilabel(app.TwoPoolTab);
            app.Pool2RefillCharacteristicTimesEditFieldLabel.HorizontalAlignment = 'right';
            app.Pool2RefillCharacteristicTimesEditFieldLabel.Position = [0 458 194 22];
            app.Pool2RefillCharacteristicTimesEditFieldLabel.Text = 'Pool 2 Refill Characteristic Time (s)';

            % Create Pool2RefillCharacteristicTimesEditField
            app.Pool2RefillCharacteristicTimesEditField = uieditfield(app.TwoPoolTab, 'numeric');
            app.Pool2RefillCharacteristicTimesEditField.Limits = [0 Inf];
            app.Pool2RefillCharacteristicTimesEditField.Position = [242 458 99 22];
            app.Pool2RefillCharacteristicTimesEditField.Value = ct2_2p;

            % Create Pool1RefillCharacteristicTimesEditFieldLabel
            app.Pool1RefillCharacteristicTimesEditFieldLabel = uilabel(app.TwoPoolTab);
            app.Pool1RefillCharacteristicTimesEditFieldLabel.HorizontalAlignment = 'right';
            app.Pool1RefillCharacteristicTimesEditFieldLabel.Position = [0 479 194 22];
            app.Pool1RefillCharacteristicTimesEditFieldLabel.Text = 'Pool 1 Refill Characteristic Time (s)';

            % Create Pool1RefillCharacteristicTimesEditField
            app.Pool1RefillCharacteristicTimesEditField = uieditfield(app.TwoPoolTab, 'numeric');
            app.Pool1RefillCharacteristicTimesEditField.Limits = [0 Inf];
            app.Pool1RefillCharacteristicTimesEditField.Position = [242 479 99 22];
            app.Pool1RefillCharacteristicTimesEditField.Value = ct1_2p;

            % Create NumberofVesiclesPool2EditFieldLabel
            app.NumberofVesiclesPool2EditFieldLabel = uilabel(app.TwoPoolTab);
            app.NumberofVesiclesPool2EditFieldLabel.HorizontalAlignment = 'right';
            app.NumberofVesiclesPool2EditFieldLabel.Position = [1 542 154 22];
            app.NumberofVesiclesPool2EditFieldLabel.Text = 'Number of Vesicles (Pool 2)';

            % Create NumberofVesiclesPool2EditField
            app.NumberofVesiclesPool2EditField = uieditfield(app.TwoPoolTab, 'numeric');
            app.NumberofVesiclesPool2EditField.Limits = [0 Inf];
            app.NumberofVesiclesPool2EditField.Position = [242 542 99 22];
            app.NumberofVesiclesPool2EditField.Value = nv2_2p;

            % Create RapidResponsePoolLabel
            app.RapidResponsePoolLabel = uilabel(app.TwoPoolTab);
            app.RapidResponsePoolLabel.Position = [350 626 122 22];
            app.RapidResponsePoolLabel.Text = 'Rapid Response Pool';

            % Create SteadyStatePoolLabel
            app.SteadyStatePoolLabel = uilabel(app.TwoPoolTab);
            app.SteadyStatePoolLabel.Position = [350 605 102 22];
            app.SteadyStatePoolLabel.Text = 'Steady State Pool';

            % Create SetParametersButton_2
            app.SetParametersButton_2 = uibutton(app.TwoPoolTab, 'push');
            app.SetParametersButton_2.ButtonPushedFcn = createCallbackFcn(app, @SetParametersButton_2Pushed, true);
            app.SetParametersButton_2.Position = [6 405 100 22];
            app.SetParametersButton_2.Text = 'Set Parameters';

            % Create SetparametersbeforegraphingLabel_2
            app.SetparametersbeforegraphingLabel_2 = uilabel(app.TwoPoolTab);
            app.SetparametersbeforegraphingLabel_2.FontWeight = 'bold';
            app.SetparametersbeforegraphingLabel_2.Position = [126 405 187 22];
            app.SetparametersbeforegraphingLabel_2.Text = 'Set parameters before graphing';

            % Create RefreshButton_2
            app.RefreshButton_2 = uibutton(app.TwoPoolTab, 'push');
            app.RefreshButton_2.ButtonPushedFcn = createCallbackFcn(app, @RefreshButton_2Pushed, true);
            app.RefreshButton_2.Position = [362 363 100 22];
            app.RefreshButton_2.Text = 'Refresh';

            % Create DependentVariableDropDown_2Label
            app.DependentVariableDropDown_2Label = uilabel(app.TwoPoolTab);
            app.DependentVariableDropDown_2Label.HorizontalAlignment = 'right';
            app.DependentVariableDropDown_2Label.Position = [0 363 111 22];
            app.DependentVariableDropDown_2Label.Text = 'Dependent Variable';

            % Create DependentVariableDropDown_2
            app.DependentVariableDropDown_2 = uidropdown(app.TwoPoolTab);
            app.DependentVariableDropDown_2.Items = {'Empty Docking Sites (Pool 1)', 'Empty Docking Sites (Pool 2)', 'Pool 1 Vesicles', 'Pool 2 Vesicles', 'Facilitated Pool 1 Vesicles', 'Facilitated Pool 2 Vesicles', 'Vesicles Released From Pool 1', 'Vesicles Released From Pool 2', 'Facilitated Vesicles Released From Pool 1', 'Facilitated Vesicles Released From Pool 2', 'Total Vesicles Released', 'Current', 'New Pool 1 Vesicle', 'New Pool 2 Vesicles', 'New Facilitated Pool 1 Vesicles', 'New Facilitated Pool 2 Vesicles'};
            app.DependentVariableDropDown_2.ItemsData = {'2', '3', '4', '5', '6', '7', '8', '9', '10', '11', '12', '13', '14', '15', '16', '17'};
            app.DependentVariableDropDown_2.ValueChangedFcn = createCallbackFcn(app, @DependentVariableDropDown_2ValueChanged, true);
            app.DependentVariableDropDown_2.Position = [126 363 224 22];
            app.DependentVariableDropDown_2.Value = '2';

            % Create UITable2
            app.UITable2 = uitable(app.TwoPoolTab);
            app.UITable2.ColumnName = {'Pulse Number'; 'Empty Docking Sites (Pool 1)'; 'Empty Docking Sites (Pool 2)'; 'Pool 1 Vesicles'; 'Pool 2 Vesicles'; 'Facilitated Pool 1 Vesicles'; 'Facilitated Pool 2 Vesicles'; 'Vesicles Released From Pool 1'; 'Vesicles Released From Pool 2'; 'Facilitated Vesicles Released From Pool 1'; 'Facilitated Vesicles Released From Pool 2'; 'Total Vesicles Released'; 'Current'; 'New Pool 1 Vesicles'; 'New Pool 2 Vesicles'; 'New Facilitated Pool 1 Vesicles'; 'New Facilitated Pool 2 Vesicles'};
            app.UITable2.RowName = {};
            app.UITable2.Position = [1 1 536 306];

            % Create RightPanel
            app.RightPanel = uipanel(app.GridLayout);
            app.RightPanel.Layout.Row = 1;
            app.RightPanel.Layout.Column = 2;

            % Create UIAxes
            app.UIAxes = uiaxes(app.RightPanel);
            title(app.UIAxes, 'Title')
            xlabel(app.UIAxes, 'X axis')
            ylabel(app.UIAxes, 'Y axis')
            app.UIAxes.Position = [1 6 657 694];

            % Show the figure after all components are created
            app.UIFigure.Visible = 'on';
            
        end
    end

    % App creation and deletion
    methods (Access = public)

        % Construct app
        function app = Synaptic_Model_Comparison_exported

            % Create UIFigure and components
            createComponents(app)

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